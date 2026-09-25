const HWSD_CN_PREPROCESSING_VERSION = v"0.2.0"
const HWSD_LAYER_BOUNDS = Float64[
    0.0 0.2
    0.2 0.4
    0.4 0.6
    0.6 0.8
    0.8 1.0
    1.0 1.5
    1.5 2.0
]
const NEURALCROP_SOIL_LAYER_BOUNDS = Float64[
    0.0 0.2
    0.2 0.5
    0.5 1.0
    1.0 2.0
    2.0 3.0
]
const _EARTH_RADIUS_M = 6_371_000.0

function _regular_spacing(coordinates::AbstractVector, label::AbstractString)
    length(coordinates) >= 2 || throw(ArgumentError("$label requires at least two coordinates"))
    differences = diff(Float64.(coordinates))
    all(>(0), differences) || all(<(0), differences) || throw(ArgumentError(
        "$label coordinates must be strictly monotonic",
    ))
    spacing = abs(first(differences))
    all(isapprox.(abs.(differences), spacing; rtol = 0, atol = max(eps(spacing) * 16, 1e-10))) ||
        throw(ArgumentError("$label coordinates must form a regular grid"))
    return spacing
end

function _coordinate_index(value::Real, coordinates::AbstractVector, spacing::Real)
    direction = coordinates[end] > coordinates[1] ? 1 : -1
    index = round(Int, (Float64(value) - coordinates[1]) / (direction * spacing)) + 1
    checkbounds(Bool, coordinates, index) || return 0
    abs(Float64(value) - coordinates[index]) <= spacing / 2 + 1e-9 || return 0
    return index
end

"""
    hwsd_tile_mapping(longitude, latitude, grid; selection=all_cells(grid))

Map the centers of one regular HWSD raster tile to compact NeuralCrop cells and
calculate each source pixel's spherical area. Returned vectors follow Julia's
`longitude × latitude` column-major order and can be passed directly to
`accumulate_soil_cn!`. Zero target indices denote pixels outside the selected
grid or over invalid canonical cells.
"""
function hwsd_tile_mapping(
    longitude::AbstractVector,
    latitude::AbstractVector,
    grid::GridIndex;
    selection::CellSelection = all_cells(grid),
)
    _validate_selection(grid, selection)
    source_lon_spacing = _regular_spacing(longitude, "source longitude")
    source_lat_spacing = _regular_spacing(latitude, "source latitude")
    target_lon_spacing = _regular_spacing(grid.longitude, "target longitude")
    target_lat_spacing = _regular_spacing(grid.latitude, "target latitude")
    lookup = zeros(Int, size(grid.cellid))
    for (local_index, compact_index) in pairs(selection.compact_indices)
        lookup[
            grid.longitude_indices[compact_index],
            grid.latitude_indices[compact_index],
        ] = local_index
    end

    pixel_count = length(longitude) * length(latitude)
    target_indices = zeros(Int32, pixel_count)
    pixel_area = zeros(Float64, pixel_count)
    position = 0
    for source_latitude in latitude, source_longitude in longitude
        position += 1
        normalized_longitude = if minimum(grid.longitude) < 0 && source_longitude > 180
            source_longitude - 360
        else
            source_longitude
        end
        longitude_index = _coordinate_index(
            normalized_longitude, grid.longitude, target_lon_spacing,
        )
        latitude_index = _coordinate_index(
            source_latitude, grid.latitude, target_lat_spacing,
        )
        if longitude_index != 0 && latitude_index != 0
            target_indices[position] = lookup[longitude_index, latitude_index]
        end
        latitude_bottom = max(-90.0, Float64(source_latitude) - source_lat_spacing / 2)
        latitude_top = min(90.0, Float64(source_latitude) + source_lat_spacing / 2)
        pixel_area[position] = _EARTH_RADIUS_M^2 * deg2rad(source_lon_spacing) *
            (sin(deg2rad(latitude_top)) - sin(deg2rad(latitude_bottom)))
    end
    return (; target_indices, pixel_area)
end

function _same_shape(arrays...)
    shape = size(first(arrays))
    all(array -> size(array) == shape, arrays) || throw(DimensionMismatch(
        "all HWSD layer attributes must have the same layer × sample shape",
    ))
    return shape
end

function _finite_or_nan(value, label)
    ismissing(value) && return NaN
    converted = Float64(value)
    isnan(converted) && return NaN
    isfinite(converted) || throw(ArgumentError("$label must be finite or missing"))
    return converted
end

"""
    hwsd_layer_stocks(organic_carbon, total_nitrogen, bulk_density, coarse_fragments)

Convert HWSD layer attributes to stocks. Inputs use the HWSD v2.x units:
organic carbon in weight percent, total nitrogen in g kg⁻¹, bulk density in
g cm⁻³, and coarse fragments in volume percent. Outputs are `(carbon,
nitrogen)` in g m⁻² for each complete source layer.
"""
function hwsd_layer_stocks(
    organic_carbon::AbstractArray,
    total_nitrogen::AbstractArray,
    bulk_density::AbstractArray,
    coarse_fragments::AbstractArray;
    layer_bounds::AbstractMatrix = HWSD_LAYER_BOUNDS,
    T::Type{<:AbstractFloat} = Float32,
)
    shape = _same_shape(
        organic_carbon, total_nitrogen, bulk_density, coarse_fragments,
    )
    size(layer_bounds) == (shape[1], 2) || throw(DimensionMismatch(
        "layer_bounds must contain one top/bottom pair per source layer",
    ))
    carbon = fill(T(NaN), shape)
    nitrogen = fill(T(NaN), shape)
    for sample in CartesianIndices(shape[2:end]), layer in axes(organic_carbon, 1)
        index = (layer, Tuple(sample)...)
        oc = _finite_or_nan(organic_carbon[index...], "organic carbon")
        tn = _finite_or_nan(total_nitrogen[index...], "total nitrogen")
        density = _finite_or_nan(bulk_density[index...], "bulk density")
        coarse = _finite_or_nan(coarse_fragments[index...], "coarse fragments")
        any(isnan, (oc, tn, density, coarse)) && continue
        0 <= oc <= 100 || throw(ArgumentError("organic carbon must be in [0, 100] wt%"))
        0 <= tn <= 1000 || throw(ArgumentError("total nitrogen must be in [0, 1000] g kg⁻¹"))
        0 < density <= 3 || throw(ArgumentError("bulk density must be in (0, 3] g cm⁻³"))
        0 <= coarse <= 100 || throw(ArgumentError("coarse fragments must be in [0, 100] vol%"))
        thickness_cm = 100 * (layer_bounds[layer, 2] - layer_bounds[layer, 1])
        thickness_cm > 0 || throw(ArgumentError("layer bounds must have positive thickness"))
        fine_fraction = 1 - coarse / 100
        carbon[index...] = T(oc * density * thickness_cm * fine_fraction * 100)
        nitrogen[index...] = T(tn * density * thickness_cm * fine_fraction * 10)
    end
    return (; carbon, nitrogen)
end

const _HWSD_ROOT_DEPTH_LAST_LAYER = (7, 5, 2, 1)

"""
    mix_hwsd_components(carbon, nitrogen, shares, root_depth)

Combine the soil components of one HWSD mapping unit without inventing deep
soil for shallow components. `root_depth` uses the HWSD classes 1–4 (deep,
moderately deep, shallow, and very shallow). Missing stocks below the final
layer represented by that class contribute zero; missing stocks inside the
represented profile remain missing. Component `SHARE` weights are never
renormalized after structural zeros are applied.
"""
function mix_hwsd_components(
    carbon::AbstractMatrix,
    nitrogen::AbstractMatrix,
    shares::AbstractVector,
    root_depth::AbstractVector;
    minimum_resolved_share::Real = 0.99,
    T::Type{<:AbstractFloat} = Float32,
)
    size(carbon) == size(nitrogen) || throw(DimensionMismatch(
        "carbon and nitrogen component stocks must have equal shape",
    ))
    size(carbon, 2) == length(shares) == length(root_depth) ||
        throw(DimensionMismatch("component metadata must match stock columns"))
    0 < minimum_resolved_share <= 1 || throw(ArgumentError(
        "minimum_resolved_share must be in (0, 1]",
    ))
    weights = T.(shares)
    all(isfinite, weights) && all(>=(zero(T)), weights) || throw(ArgumentError(
        "component shares must be finite and nonnegative",
    ))
    total_share = sum(weights)
    total_share > zero(T) || throw(ArgumentError("component shares must sum to a positive value"))
    weights ./= total_share

    mixed_carbon = fill(T(NaN), size(carbon, 1))
    mixed_nitrogen = fill(T(NaN), size(nitrogen, 1))
    uncertain = falses(size(carbon, 1))
    for layer in axes(carbon, 1)
        carbon_total = zero(T)
        nitrogen_total = zero(T)
        resolved_share = zero(T)
        for component in axes(carbon, 2)
            carbon_value = carbon[layer, component]
            nitrogen_value = nitrogen[layer, component]
            root_class = Int(root_depth[component])
            structural_zero = 1 <= root_class <= length(_HWSD_ROOT_DEPTH_LAST_LAYER) &&
                layer > _HWSD_ROOT_DEPTH_LAST_LAYER[root_class]
            if isfinite(carbon_value) && isfinite(nitrogen_value)
                carbon_total += weights[component] * T(carbon_value)
                nitrogen_total += weights[component] * T(nitrogen_value)
                resolved_share += weights[component]
            elseif structural_zero
                resolved_share += weights[component]
                uncertain[layer] = true
            end
        end
        if resolved_share >= T(minimum_resolved_share)
            mixed_carbon[layer] = carbon_total
            mixed_nitrogen[layer] = nitrogen_total
            uncertain[layer] |= resolved_share < one(T)
        end
    end
    return (; carbon = mixed_carbon, nitrogen = mixed_nitrogen, uncertain)
end

function _validate_layer_bounds(bounds::AbstractMatrix, label::AbstractString)
    size(bounds, 2) == 2 || throw(DimensionMismatch("$label must have two columns"))
    all(bounds[:, 2] .> bounds[:, 1]) || throw(ArgumentError(
        "$label layers must have positive thickness",
    ))
    all(bounds[2:end, 1] .== bounds[1:end-1, 2]) || throw(ArgumentError(
        "$label layers must be contiguous",
    ))
    return nothing
end

"""
    remap_hwsd_layers(stocks; deep_rule=:extend_deepest_density)

Map complete-layer stocks from the seven HWSD layers to the five NeuralCrop
layers by exact vertical overlap. HWSD stops at 2 m. The default rule extends
the 1.5–2 m stock density through 2–3 m and flags that target layer as
uncertain. Use `deep_rule=:missing` to leave unsupported depth missing.
"""
function remap_hwsd_layers(
    stocks::AbstractMatrix;
    source_bounds::AbstractMatrix = HWSD_LAYER_BOUNDS,
    target_bounds::AbstractMatrix = NEURALCROP_SOIL_LAYER_BOUNDS,
    deep_rule::Symbol = :extend_deepest_density,
    T::Type{<:AbstractFloat} = eltype(stocks) <: AbstractFloat ? eltype(stocks) : Float32,
)
    _validate_layer_bounds(source_bounds, "source_bounds")
    _validate_layer_bounds(target_bounds, "target_bounds")
    size(stocks, 1) == size(source_bounds, 1) || throw(DimensionMismatch(
        "stock rows must match source layers",
    ))
    deep_rule in (:extend_deepest_density, :missing) || throw(ArgumentError(
        "deep_rule must be :extend_deepest_density or :missing",
    ))
    samples = size(stocks, 2)
    remapped = fill(T(NaN), size(target_bounds, 1), samples)
    uncertain = falses(size(remapped))
    source_bottom = source_bounds[end, 2]
    for sample in 1:samples, target in axes(target_bounds, 1)
        top, bottom = target_bounds[target, :]
        total = zero(T)
        valid = true
        for source in axes(source_bounds, 1)
            overlap = max(0.0, min(bottom, source_bounds[source, 2]) -
                max(top, source_bounds[source, 1]))
            overlap == 0 && continue
            value = stocks[source, sample]
            if !isfinite(value)
                valid = false
                break
            end
            source_thickness = source_bounds[source, 2] - source_bounds[source, 1]
            total += T(value * overlap / source_thickness)
        end
        unsupported = max(0.0, bottom - max(top, source_bottom))
        if unsupported > 0
            uncertain[target, sample] = true
            if deep_rule === :missing || !isfinite(stocks[end, sample])
                valid = false
            else
                deepest_thickness = source_bounds[end, 2] - source_bounds[end, 1]
                total += T(stocks[end, sample] * unsupported / deepest_thickness)
            end
        end
        valid && (remapped[target, sample] = total)
    end
    return (; values = remapped, uncertain)
end

function init_soil_cn_aggregator(
    ::Type{T}, layers::Integer, cells::Integer,
) where {T <: AbstractFloat}
    layers > 0 || throw(ArgumentError("layers must be positive"))
    cells > 0 || throw(ArgumentError("cells must be positive"))
    matrix() = zeros(T, layers, cells)
    return SoilCNAggregator(
        matrix(), matrix(), matrix(), matrix(), matrix(), zeros(T, cells),
    )
end

init_soil_cn_aggregator(layers::Integer, cells::Integer) =
    init_soil_cn_aggregator(Float64, layers, cells)

"""
    accumulate_soil_cn!(accumulator, carbon, nitrogen, target_indices, pixel_area)

Add one source-raster tile. `carbon` and `nitrogen` are layer × pixel stocks in
g m⁻²; each pixel maps to a one-based compact target index, or zero outside the
target grid. This incremental interface keeps preprocessing memory bounded by
one HWSD tile.
"""
function accumulate_soil_cn!(
    accumulator::SoilCNAggregator{T},
    carbon::AbstractMatrix,
    nitrogen::AbstractMatrix,
    target_indices::AbstractVector{<:Integer},
    pixel_area::AbstractVector;
    uncertain::AbstractMatrix{Bool} = falses(size(carbon)),
) where {T <: AbstractFloat}
    size(carbon) == size(nitrogen) == size(uncertain) || throw(DimensionMismatch(
        "carbon, nitrogen, and uncertainty arrays must have equal shape",
    ))
    size(carbon, 1) == size(accumulator.carbon_sum, 1) || throw(DimensionMismatch(
        "source target-layer count does not match accumulator",
    ))
    size(carbon, 2) == length(target_indices) == length(pixel_area) ||
        throw(DimensionMismatch("pixel arrays must have equal length"))
    for pixel in eachindex(target_indices)
        target = Int(target_indices[pixel])
        target == 0 && continue
        checkbounds(Bool, accumulator.target_area, target) || throw(BoundsError(
            accumulator.target_area, target,
        ))
        area = T(pixel_area[pixel])
        isfinite(area) && area > 0 || throw(ArgumentError("pixel areas must be positive and finite"))
        accumulator.target_area[target] += area
        for layer in axes(carbon, 1)
            carbon_value = carbon[layer, pixel]
            if isfinite(carbon_value)
                accumulator.carbon_sum[layer, target] += T(carbon_value) * area
                accumulator.carbon_area[layer, target] += area
            end
            nitrogen_value = nitrogen[layer, pixel]
            if isfinite(nitrogen_value)
                accumulator.nitrogen_sum[layer, target] += T(nitrogen_value) * area
                accumulator.nitrogen_area[layer, target] += area
            end
            uncertain[layer, pixel] &&
                (accumulator.uncertain_area[layer, target] += area)
        end
    end
    return accumulator
end

function finish_soil_cn(
    accumulator::SoilCNAggregator{T},
    selection::CellSelection;
    layer_bounds::AbstractMatrix = NEURALCROP_SOIL_LAYER_BOUNDS,
    minimum_coverage::Real = 0.99,
    provenance::NamedTuple,
) where {T <: AbstractFloat}
    0 < minimum_coverage <= 1 || throw(ArgumentError("minimum_coverage must be in (0, 1]"))
    cells = length(selection.cell_ids)
    cells == length(accumulator.target_area) || throw(DimensionMismatch(
        "selection does not match accumulator cells",
    ))
    layers = size(accumulator.carbon_sum, 1)
    size(layer_bounds) == (layers, 2) || throw(DimensionMismatch(
        "layer bounds do not match accumulator layers",
    ))
    carbon = fill(T(NaN), layers, cells)
    nitrogen = fill(T(NaN), layers, cells)
    coverage = zeros(T, layers, cells)
    uncertain = falses(layers, cells)
    for cell in 1:cells, layer in 1:layers
        target_area = accumulator.target_area[cell]
        target_area == 0 && continue
        carbon_coverage = accumulator.carbon_area[layer, cell] / target_area
        nitrogen_coverage = accumulator.nitrogen_area[layer, cell] / target_area
        joint_coverage = min(carbon_coverage, nitrogen_coverage)
        coverage[layer, cell] = joint_coverage
        if joint_coverage >= minimum_coverage
            carbon[layer, cell] = accumulator.carbon_sum[layer, cell] /
                accumulator.carbon_area[layer, cell]
            nitrogen[layer, cell] = accumulator.nitrogen_sum[layer, cell] /
                accumulator.nitrogen_area[layer, cell]
        end
        uncertain[layer, cell] = joint_coverage < one(T) ||
            accumulator.uncertain_area[layer, cell] > zero(T)
    end
    return SoilCNTargets(
        selection, T.(layer_bounds), carbon, nitrogen, coverage, uncertain,
        provenance,
    )
end

"""Return layer-wise source totals represented by an area accumulator."""
function soil_cn_conservation(accumulator::SoilCNAggregator)
    return (
        carbon_g = vec(sum(accumulator.carbon_sum; dims = 2)),
        nitrogen_g = vec(sum(accumulator.nitrogen_sum; dims = 2)),
        carbon_valid_area_m2 = vec(sum(accumulator.carbon_area; dims = 2)),
        nitrogen_valid_area_m2 = vec(sum(accumulator.nitrogen_area; dims = 2)),
    )
end

"""Run the numerical HWSD C/N conversion for one or more streamed raster tiles."""
function preprocess_hwsd_cn(
    organic_carbon::AbstractMatrix,
    total_nitrogen::AbstractMatrix,
    bulk_density::AbstractMatrix,
    coarse_fragments::AbstractMatrix,
    target_indices::AbstractVector{<:Integer},
    pixel_area::AbstractVector,
    selection::CellSelection;
    source_version::AbstractString = "HWSD v2.01",
    deep_rule::Symbol = :extend_deepest_density,
    minimum_coverage::Real = 0.99,
    T::Type{<:AbstractFloat} = Float32,
)
    stocks = hwsd_layer_stocks(
        organic_carbon, total_nitrogen, bulk_density, coarse_fragments; T,
    )
    carbon = remap_hwsd_layers(stocks.carbon; deep_rule, T)
    nitrogen = remap_hwsd_layers(stocks.nitrogen; deep_rule, T)
    aggregator = init_soil_cn_aggregator(T, size(NEURALCROP_SOIL_LAYER_BOUNDS, 1), length(selection.cell_ids))
    accumulate_soil_cn!(
        aggregator, carbon.values, nitrogen.values, target_indices, pixel_area;
        uncertain = carbon.uncertain .| nitrogen.uncertain,
    )
    provenance = (
        schema_version = DATA_SCHEMA_VERSION,
        preprocessing_version = HWSD_CN_PREPROCESSING_VERSION,
        source_version = String(source_version),
        source_units = (
            organic_carbon = "wt%", total_nitrogen = "g/kg",
            bulk_density = "g/cm3", coarse_fragments = "vol%",
        ),
        output_units = (carbon = "gC/m2", nitrogen = "gN/m2"),
        deep_rule,
        conservation = soil_cn_conservation(aggregator),
    )
    return finish_soil_cn(
        aggregator, selection; minimum_coverage, provenance,
    )
end

function write_soil_cn_targets(path::AbstractString, targets::SoilCNTargets)
    isfile(path) && rm(path; force = true)
    NCDataset(path, "c") do dataset
        layers, cells = size(targets.soil_organic_carbon)
        defDim(dataset, "layer", layers)
        defDim(dataset, "bounds", 2)
        defDim(dataset, "cell", cells)
        defVar(dataset, "cell_id", Int32, ("cell",))[:] = targets.selection.cell_ids
        defVar(dataset, "layer_bounds", eltype(targets.layer_bounds), ("layer", "bounds"))[:, :] =
            targets.layer_bounds
        carbon = defVar(dataset, "soil_organic_carbon", eltype(targets.soil_organic_carbon), ("layer", "cell"))
        nitrogen = defVar(dataset, "total_nitrogen", eltype(targets.total_nitrogen), ("layer", "cell"))
        coverage = defVar(dataset, "coverage", eltype(targets.coverage), ("layer", "cell"))
        uncertain = defVar(dataset, "uncertain", Int8, ("layer", "cell"))
        carbon[:, :] = targets.soil_organic_carbon
        nitrogen[:, :] = targets.total_nitrogen
        coverage[:, :] = targets.coverage
        uncertain[:, :] = Int8.(targets.uncertain)
        carbon.attrib["units"] = "gC m-2"
        nitrogen.attrib["units"] = "gN m-2"
        coverage.attrib["units"] = "1"
        dataset.attrib["schema_version"] = string(targets.provenance.schema_version)
        dataset.attrib["preprocessing_version"] = string(targets.provenance.preprocessing_version)
        dataset.attrib["source_version"] = targets.provenance.source_version
        dataset.attrib["deep_rule"] = string(targets.provenance.deep_rule)
        for name in (
            :fill_policy, :donor_longitude, :donor_latitude,
            :fill_distance_km, :original_minimum_coverage,
        )
            haskey(targets.provenance, name) || continue
            dataset.attrib[string(name)] = targets.provenance[name]
        end
    end
    return String(path)
end

function read_soil_cn_targets(path::AbstractString; T::Type{<:AbstractFloat} = Float32)
    return NCDataset(path, "r") do dataset
        cell_ids = Int32.(dataset["cell_id"][:])
        selection = CellSelection(eachindex(cell_ids), cell_ids)
        provenance = (
            schema_version = VersionNumber(dataset.attrib["schema_version"]),
            preprocessing_version = VersionNumber(dataset.attrib["preprocessing_version"]),
            source_version = String(dataset.attrib["source_version"]),
            deep_rule = Symbol(dataset.attrib["deep_rule"]),
        )
        for (name, convert) in (
            (:fill_policy, String),
            (:donor_longitude, Float64),
            (:donor_latitude, Float64),
            (:fill_distance_km, Float64),
            (:original_minimum_coverage, Float64),
        )
            attribute = string(name)
            haskey(dataset.attrib, attribute) || continue
            provenance = merge(provenance, NamedTuple{(name,)}((convert(dataset.attrib[attribute]),)))
        end
        return SoilCNTargets(
            selection,
            T.(dataset["layer_bounds"][:, :]),
            T.(dataset["soil_organic_carbon"][:, :]),
            T.(dataset["total_nitrogen"][:, :]),
            T.(dataset["coverage"][:, :]),
            BitMatrix(dataset["uncertain"][:, :] .!= 0),
            provenance,
        )
    end
end

"""Write one CFT/water-system-calibrated soil-pool allocation for later reconstruction."""
function write_soil_pool_allocation(path::AbstractString, allocation::SoilPoolAllocation)
    isfile(path) && rm(path; force = true)
    NCDataset(path, "c") do dataset
        layers, cells = size(allocation.fast_carbon_fraction)
        defDim(dataset, "layer", layers)
        defDim(dataset, "cell", cells)
        defDim(dataset, "patch", 1)
        defVar(dataset, "cell_id", Int32, ("cell",))[:] = allocation.selection.cell_ids
        defVar(dataset, "cft_id", Int32, ("patch",))[:] = Int32[allocation.cft_id]
        defVar(dataset, "irrigated", Int8, ("patch",))[:] = Int8[allocation.irrigated]
        for name in (
            :fast_carbon_fraction,
            :fast_nitrogen_fraction,
            :c_shift_fast,
            :c_shift_slow,
        )
            values = getproperty(allocation, name)
            defVar(dataset, String(name), eltype(values), ("layer", "cell", "patch"))[:, :, 1] = values
        end
        dataset.attrib["schema_version"] = string(DATA_SCHEMA_VERSION)
        for (name, value) in pairs(allocation.provenance)
            value isa Union{AbstractString, Number, Symbol} || continue
            dataset.attrib[string(name)] = string(value)
        end
    end
    return String(path)
end

"""Read one CFT/water-system-calibrated soil-pool allocation."""
function read_soil_pool_allocation(path::AbstractString; T::Type{<:AbstractFloat} = Float32)
    return NCDataset(path, "r") do dataset
        cell_ids = Int32.(dataset["cell_id"][:])
        selection = CellSelection(eachindex(cell_ids), cell_ids)
        provenance = (schema_version = VersionNumber(dataset.attrib["schema_version"]),)
        for name in keys(dataset.attrib)
            name == "schema_version" && continue
            provenance = merge(provenance, NamedTuple{(Symbol(name),)}((dataset.attrib[name],)))
        end
        # `pft_id` was written by pre-CFT allocation files. Accept it only at
        # this serialized-data boundary and normalize immediately to `cft_id`.
        cft_variable = haskey(dataset, "cft_id") ? "cft_id" :
            (haskey(dataset, "pft_id") ? "pft_id" : nothing)
        has_patch = !isnothing(cft_variable) && haskey(dataset, "irrigated")
        has_patch && length(dataset[cft_variable]) != 1 && error(
            "allocation file contains multiple patches; select one CFT/water-system allocation first",
        )
        cft_id = has_patch ? Int(dataset[cft_variable][1]) : 1
        irrigated = has_patch ? Bool(dataset["irrigated"][1]) : false
        values(name) = ndims(dataset[name]) == 3 ? dataset[name][:, :, 1] : dataset[name][:, :]
        return SoilPoolAllocation(
            selection,
            T.(values("fast_carbon_fraction")),
            T.(values("fast_nitrogen_fraction")),
            T.(values("c_shift_fast")),
            T.(values("c_shift_slow"));
            cft_id, irrigated,
            provenance,
        )
    end
end
