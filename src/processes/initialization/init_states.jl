"""
root_distribution(beta_root)

Compute normalized root fractions across the default five soil layers from
the LPJmL-style exponential root profile parameter `beta_root`.
"""
function root_distribution(beta_root::AbstractFloat)
    T = typeof(beta_root)
    layerbound = T[200.0, 500.0, 1000.0, 2000.0, 3000.0]

    BOTTOMLAYER = length(layerbound)
    totalroots = one(T) - beta_root^(layerbound[BOTTOMLAYER] / T(10))
    rootdist = zeros(T, BOTTOMLAYER)
    rootdist[1] = (one(T) - beta_root^(layerbound[1] / T(10))) / totalroots
    for l in 2:BOTTOMLAYER
        rootdist[l] = (
            beta_root^(layerbound[l-1] / T(10)) -
            beta_root^(layerbound[l] / T(10))
        ) / totalroots
    end

    return rootdist
end

@kernel inbounds = true function initialize_mineral_nitrogen_kernel!(
    nitrate, ammonium, slow, nitrate_restart, ammonium_restart,
    slow_fraction, restart::Bool,
)
    layer, cell = @index(Global, NTuple)
    if restart
        nitrate[layer, cell] = nitrate_restart[layer, cell]
        ammonium[layer, cell] = ammonium_restart[layer, cell]
    else
        value = slow[layer, cell] * slow_fraction
        nitrate[layer, cell] = value
        ammonium[layer, cell] = value
    end
end

@kernel inbounds = true function initialize_crop_phenology_kernel!(
    root_distribution_state, phu, winter_type, root_distribution_values,
)
    layer, cell = @index(Global, NTuple)
    root_distribution_state[layer, cell] = root_distribution_values[layer]
    if layer == 1 && phu[cell] < zero(eltype(phu))
        phu[cell] = -phu[cell]
        winter_type[cell] = true
    end
end

@kernel inbounds = true function initialize_c_shift_kernel!(
    shift_fast, shift_slow, restart_fast, restart_slow,
    top_fraction, lower_fraction, restart::Bool,
)
    layer, cell = @index(Global, NTuple)
    if restart
        shift_fast[layer, cell] = restart_fast[layer, cell]
        shift_slow[layer, cell] = restart_slow[layer, cell]
    else
        value = layer == 1 ? top_fraction : lower_fraction
        shift_fast[layer, cell] = value
        shift_slow[layer, cell] = value
    end
end

"""
initialize_soil_mineral_nitrogen!(soil, u0, strategy)

Initialize soil NO₃ and NH₄ using either restart values or the default rule.
The default initializes each mineral pool in each layer to one percent of that
layer's slow organic-N pool.
"""
function initialize_soil_mineral_nitrogen!(soil::Soil,
                                           u0::NamedTuple,
                                           strategy::Symbol)
    restart = strategy === :restart
    if strategy === :restart
        if !hasproperty(u0, :soil_NO3) || !hasproperty(u0, :soil_NH4)
            throw(ArgumentError(
                ":restart requires soil_NO3 and soil_NH4; construct inputs " *
                "with load_mineral_nitrogen_restart=true",
            ))
        end
    elseif strategy === :from_slow_organic_nitrogen
    else
        throw(ArgumentError(
            "unknown mineral nitrogen initialization strategy: $strategy; " *
            "use :restart or :from_slow_organic_nitrogen",
        ))
    end
    slow_n_fraction = convert(eltype(soil.nitrogen.slow), 0.01)
    nitrate_restart = restart ? u0.soil_NO3 : soil.nitrogen.nitrate
    ammonium_restart = restart ? u0.soil_NH4 : soil.nitrogen.ammonium
    launch_2D!(
        initialize_mineral_nitrogen_kernel!, soil.nitrogen.nitrate,
        soil.nitrogen.ammonium, soil.nitrogen.slow,
        nitrate_restart, ammonium_restart, slow_n_fraction, restart,
    )
    return nothing
end


"""
init_states!(CFT, InitialData, cell_size, device;
             lpjmlparams=lpjmlparams,
             mineral_nitrogen_initialization=:from_slow_organic_nitrogen)

Initialize and populate all runtime state structs from static parameters and
input data for one simulation domain.
Returns `(climbuf, crop, pet, soil, managed_land, dailyWeather, output)`.
Crop storage is separated by lifetime into `crop.state`, `crop.fluxes`,
`crop.auxiliary`, and `crop.workspace`.
"""
function init_states!(CFT::CFTParameters,
                       InitialData::NamedTuple,
                       cell_size::Int,
                       device;
                       T::Type{<:AbstractFloat} = Float32,
                       lpjmlparams::LPJmLParams = lpjmlparams,
                       mineral_nitrogen_initialization::Symbol = :from_slow_organic_nitrogen,
                       c_shift_initialization::Symbol = :default
)

    @unpack residue_frac = lpjmlparams
    @unpack k_litter10, beta_root = CFT

    @unpack latitude, soilparams, ModelState = InitialData

    phu = ModelState.crop.phu
    sdate = ModelState.crop.sdate
    manure = ModelState.crop.manure
    fertilizer = ModelState.crop.fertilizer
    residuefrac = ModelState.crop.residuefrac
    u0 = ModelState.u0

    to_float(values) = device(T.(values))
    to_integer(values) = device(Int32.(values))

    dailyWeather = init_weather(T, cell_size, device)
    climbuf = init_climbuf(T, cell_size, device)
    crop = init_crop(T, cell_size, device)
    managed_land = init_managed_land(T, cell_size, device)
    crop.auxiliary.phenology.phu = to_float(phu)
    rootdist = root_distribution(T(beta_root))
    launch_2D!(
        initialize_crop_phenology_kernel!, crop.auxiliary.root.distribution,
        crop.auxiliary.phenology.phu, crop.auxiliary.phenology.winter_type,
        device(rootdist),
    )

    crop.auxiliary.calendar.sowing_date = to_integer(sdate)
    crop.auxiliary.calendar.prescribed_sowing_date = to_integer(sdate)
    managed_land.manure = to_float(manure)
    managed_land.fertilizer = to_float(fertilizer)
    managed_land.residue_fraction = to_float(residuefrac)
    managed_land.latitude = to_float(latitude)
    pet = init_pet(T, cell_size, device)
    soil = init_soil(T, cell_size, T.(soilparams.soildepth), device)
    soil.carbon.litter = to_float(u0.litc)
    soil.carbon.fast = to_float(u0.fastc)
    soil.carbon.slow = to_float(u0.slowc)
    soil.nitrogen.litter = to_float(u0.litn)
    soil.nitrogen.fast = to_float(u0.fastn)
    soil.nitrogen.slow = to_float(u0.slown)
    soil.water.storage = to_float(u0.swc)
    mineral_u0 = mineral_nitrogen_initialization === :restart ? merge(u0, (
        soil_NO3 = to_float(u0.soil_NO3),
        soil_NH4 = to_float(u0.soil_NH4),
    )) : u0
    initialize_soil_mineral_nitrogen!(
        soil,
        mineral_u0,
        mineral_nitrogen_initialization,
    )
    soil.water.saturation_fraction = to_float(soilparams.w_sat)
    soil.properties.ph = to_float(soilparams.ph)
    soil.properties.sand_fraction = to_float(soilparams.sand)
    soil.properties.clay_fraction = to_float(soilparams.clay)
    soil.thermal.diffusivity_0 = to_float(soilparams.tdiff_0)
    soil.thermal.diffusivity_15 = to_float(soilparams.tdiff_15)

    soil.management.tillage_fraction = device(T[
        1 - residue_frac 0 0
        residue_frac 1 0
        0 0 1
    ])
    c_shift_state = c_shift_initialization === :restart ? merge(ModelState, (
        c_shift_fast = to_float(ModelState.c_shift_fast),
        c_shift_slow = to_float(ModelState.c_shift_slow),
    )) : ModelState
    initialize_soil_c_shift!(soil, c_shift_state, c_shift_initialization)
    days_per_year = T(365)
    soil.carbon.litter_response = device(
        T[k_litter10.leaf, k_litter10.leaf, k_litter10.root] ./ days_per_year,
    )

    soil.nitrogen.litter_response = device(
        T[k_litter10.leaf, k_litter10.leaf, k_litter10.root] ./ days_per_year,
    )

    output = init_output(T, cell_size, device)

    synchronize_backend!(crop.auxiliary.phenology.phu)

    return climbuf, crop, pet, soil, managed_land, dailyWeather, output
end

"""
    initialize_soil_c_shift!(soil, model_state, strategy)

Initialize the normalized vertical distribution used to route decomposed
litter into fast and slow soil pools.

`:default` uses 0.55 in the top layer and 0.45 distributed uniformly over all
remaining layers. `:restart` restores the fast and slow distributions supplied
in `model_state`.
"""
function initialize_soil_c_shift!(soil::Soil,
                                  model_state::NamedTuple,
                                  strategy::Symbol)
    if strategy === :default
        layers = size(soil.decomposition.shift_fast, 1)
        layers > 1 || throw(ArgumentError("default c_shift initialization requires at least two soil layers"))
        T = eltype(soil.decomposition.shift_fast)
        lower_layer_fraction = T(0.45) / T(layers - 1)

        restart_fast = soil.decomposition.shift_fast
        restart_slow = soil.decomposition.shift_slow
        restart = false
    elseif strategy === :restart
        hasproperty(model_state, :c_shift_fast) ||
            throw(ArgumentError("c_shift_initialization=:restart requires ModelState.c_shift_fast"))
        hasproperty(model_state, :c_shift_slow) ||
            throw(ArgumentError("c_shift_initialization=:restart requires ModelState.c_shift_slow"))
        size(model_state.c_shift_fast) == size(soil.decomposition.shift_fast) ||
            throw(DimensionMismatch("c_shift_fast must match the soil layer-by-cell shape"))
        size(model_state.c_shift_slow) == size(soil.decomposition.shift_slow) ||
            throw(DimensionMismatch("c_shift_slow must match the soil layer-by-cell shape"))
        T = eltype(soil.decomposition.shift_fast)
        lower_layer_fraction = zero(T)
        restart_fast = model_state.c_shift_fast
        restart_slow = model_state.c_shift_slow
        restart = true
    else
        throw(ArgumentError("unknown c_shift initialization strategy: $strategy"))
    end

    launch_2D!(
        initialize_c_shift_kernel!, soil.decomposition.shift_fast,
        soil.decomposition.shift_slow, restart_fast, restart_slow,
        T(0.55), lower_layer_fraction, restart,
    )

    return nothing
end
