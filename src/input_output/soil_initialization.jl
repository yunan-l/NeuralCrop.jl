const _INITIAL_SOIL_STATE_FIELDS = (:swc, :litc, :fastc, :slowc, :litn, :fastn, :slown)

"""Calculate initial liquid water storage at field capacity in mm."""
function field_capacity_water(
    soil,
    soil_organic_carbon::AbstractMatrix;
    mineral_density::Real = 2700,
    T::Type{<:AbstractFloat} = Float32,
)
    size(soil_organic_carbon) == size(soil.saturation) || throw(DimensionMismatch(
        "soil organic carbon must match the soil layer × cell shape",
    ))
    sand = reshape(T.(soil.sand), 1, :)
    clay = reshape(T.(soil.clay), 1, :)
    depth = reshape(T.(soil.layer_depth), :, 1)
    saturation_reference = T.(soil.saturation)
    carbon = T.(soil_organic_carbon)
    organic_matter = clamp.(
        T(2) .* carbon ./
            ((one(T) .- saturation_reference) .* T(mineral_density) .* depth) .* T(100),
        zero(T), T(8),
    )
    ws33t = T(0.278) .* sand .+ T(0.034) .* clay .+
        T(0.022) .* organic_matter .- T(0.018) .* sand .* organic_matter .-
        T(0.027) .* clay .* organic_matter .-
        T(0.584) .* sand .* clay .+ T(0.078)
    ws33 = ws33t .+ (T(0.636) .* ws33t .- T(0.107))
    wfct = -T(0.251) .* sand .+ T(0.195) .* clay .+
        T(0.011) .* organic_matter .+ T(0.006) .* sand .* organic_matter .-
        T(0.027) .* clay .* organic_matter .+
        T(0.452) .* sand .* clay .+ T(0.299)
    field = wfct .+ ((T(1.283) .* wfct) .^ 2 .- T(0.374) .* wfct .- T(0.015))
    saturation = field .+ ws33 .- T(0.097) .* sand .+ T(0.043)
    field = ifelse.(saturation .- field .< T(0.05), saturation .- T(0.05), field)
    return field .* depth
end

"""
    soil_initial_state(targets, soil; fast_fraction=0.4, allocation=nothing)

Build NeuralCrop's seven initial state arrays from externally supplied, aligned
soil C/N targets and static soil properties. This is model initialization
policy, not input-data decoding. The default 40:60 split is replaced by a
supplied neutral pool-allocation contract when available.
"""
function soil_initial_state(
    targets,
    soil;
    fast_fraction::Real = 0.4,
    allocation = nothing,
    mineral_fraction_of_slow::Real = 0.01,
    T::Type{<:AbstractFloat} = Float32,
)
    targets.selection.cell_ids == soil.selection.cell_ids || throw(ArgumentError(
        "soil targets and soil properties must use the same compact cell selection",
    ))
    size(targets.soil_organic_carbon) == size(soil.saturation) || throw(DimensionMismatch(
        "soil targets and soil properties must have matching layers",
    ))
    mineral_fraction_of_slow >= 0 || throw(ArgumentError(
        "mineral_fraction_of_slow must be nonnegative",
    ))
    carbon = T.(targets.soil_organic_carbon)
    nitrogen = T.(targets.total_nitrogen)
    all(isfinite, carbon) || throw(ArgumentError("soil SOC targets contain missing values"))
    all(isfinite, nitrogen) || throw(ArgumentError("soil total-N targets contain missing values"))
    all(>=(zero(T)), carbon) || throw(ArgumentError("soil SOC targets must be nonnegative"))
    all(>=(zero(T)), nitrogen) || throw(ArgumentError("soil total-N targets must be nonnegative"))

    fast_carbon, fast_nitrogen, c_shift = if allocation === nothing
        0 <= fast_fraction <= 1 || throw(ArgumentError("fast_fraction must be in [0, 1]"))
        fast = fill(T(fast_fraction), size(carbon))
        (fast, copy(fast), nothing)
    else
        allocation.selection.cell_ids == targets.selection.cell_ids || throw(ArgumentError(
            "pool allocation and HWSD targets must use the same compact cell selection",
        ))
        for name in (:fast_carbon_fraction, :fast_nitrogen_fraction, :c_shift_fast, :c_shift_slow)
            values = getproperty(allocation, name)
            size(values) == size(carbon) || throw(DimensionMismatch(
                "pool allocation.$name must match the HWSD layer × cell shape",
            ))
            all(isfinite, values) || throw(ArgumentError("pool allocation.$name must be finite"))
        end
        all((zero(eltype(allocation.fast_carbon_fraction)) .<= allocation.fast_carbon_fraction) .&
            (allocation.fast_carbon_fraction .<= one(eltype(allocation.fast_carbon_fraction)))) ||
            throw(ArgumentError("fast carbon fractions must be in [0, 1]"))
        all((zero(eltype(allocation.fast_nitrogen_fraction)) .<= allocation.fast_nitrogen_fraction) .&
            (allocation.fast_nitrogen_fraction .<= one(eltype(allocation.fast_nitrogen_fraction)))) ||
            throw(ArgumentError("fast nitrogen fractions must be in [0, 1]"))
        for name in (:c_shift_fast, :c_shift_slow)
            values = getproperty(allocation, name)
            all(values .>= zero(eltype(values))) || throw(ArgumentError(
                "pool allocation.$name must be nonnegative",
            ))
            all(isapprox.(vec(sum(values; dims = 1)), one(eltype(values)); atol = sqrt(eps(eltype(values))))) ||
                throw(ArgumentError("pool allocation.$name must sum to one by cell"))
        end
        (
            T.(allocation.fast_carbon_fraction),
            T.(allocation.fast_nitrogen_fraction),
            (fast = T.(allocation.c_shift_fast), slow = T.(allocation.c_shift_slow)),
        )
    end
    slow_carbon = one(T) .- fast_carbon
    slow_nitrogen = one(T) .- fast_nitrogen
    mineral = T(mineral_fraction_of_slow)
    organic_nitrogen = nitrogen ./ (one(T) .+ T(2) * mineral .* slow_nitrogen)
    cells = length(targets.selection.cell_ids)
    state = (
        swc = field_capacity_water(soil, carbon; T),
        litc = zeros(T, 3, cells),
        fastc = fast_carbon .* carbon,
        slowc = slow_carbon .* carbon,
        litn = zeros(T, 3, cells),
        fastn = fast_nitrogen .* organic_nitrogen,
        slown = slow_nitrogen .* organic_nitrogen,
    )
    return c_shift === nothing ? state : merge(state, (
        c_shift_fast = c_shift.fast,
        c_shift_slow = c_shift.slow,
    ))
end

"""Build backend-neutral initial data consumed by `initialize_simulation`."""
function model_initial_data(grid, soil, crop::NamedTuple, initial_state::NamedTuple)
    selection = soil.selection
    cells = length(selection.cell_ids)
    for name in (:sdate, :phu, :manure, :fertilizer, :residuefrac)
        hasproperty(crop, name) || throw(ArgumentError("crop inputs are missing $name"))
        length(getproperty(crop, name)) == cells || throw(DimensionMismatch(
            "crop.$name must contain one value per selected cell",
        ))
    end
    for name in _INITIAL_SOIL_STATE_FIELDS
        hasproperty(initial_state, name) || throw(ArgumentError("initial_state is missing $name"))
        values = getproperty(initial_state, name)
        ndims(values) == 2 || throw(DimensionMismatch("initial_state.$name must be layer × cell"))
        size(values, 2) == cells || throw(DimensionMismatch(
            "initial_state.$name must use the soil cell selection",
        ))
    end
    has_fast_shift = hasproperty(initial_state, :c_shift_fast)
    has_slow_shift = hasproperty(initial_state, :c_shift_slow)
    has_fast_shift == has_slow_shift || throw(ArgumentError(
        "initial_state must provide both c_shift_fast and c_shift_slow",
    ))
    if has_fast_shift
        for name in (:c_shift_fast, :c_shift_slow)
            size(getproperty(initial_state, name)) == size(soil.saturation) || throw(
                DimensionMismatch("initial_state.$name must use the soil layer × cell shape"),
            )
        end
    end
    compact = selection.compact_indices
    data = (
        coords = copy(selection.cell_ids),
        latitude = grid.latitude[grid.latitude_indices[compact]],
        crop,
        soilparam = (
            soilcode = soil.soilcode,
            soilph = soil.ph,
            w_sat = soil.saturation,
            sand = soil.sand,
            silt = soil.silt,
            clay = soil.clay,
            tdiff_0 = soil.diffusivity_dry,
            tdiff_15 = soil.diffusivity_15,
            soildepth = soil.layer_depth,
        ),
        initial_state,
        backend_neutral = true,
    )
    return has_fast_shift ? merge(data, (
        c_shift_fast = copy(initial_state.c_shift_fast),
        c_shift_slow = copy(initial_state.c_shift_slow),
    )) : data
end
