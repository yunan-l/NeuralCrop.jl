const DEFAULT_SOIL_LOOKUP_VERSION = v"1.0.0"

"""Return the versioned 14-class soil lookup currently used by NeuralCrop."""
function default_soil_lookup(::Type{T} = Float32) where {T <: AbstractFloat}
    return SoilLookup(
        DEFAULT_SOIL_LOOKUP_VERSION,
        T[0.22, 0.06, 0.52, 0.32, 0.10, 0.58, 0.43, 0.17, 0.58, 0.10, 0.82, 0.92, 0.24, 0.99],
        T[0.20, 0.47, 0.06, 0.34, 0.56, 0.15, 0.39, 0.70, 0.32, 0.60, 0.12, 0.05, 0.28, 0.00],
        T[0.58, 0.47, 0.42, 0.34, 0.34, 0.27, 0.18, 0.13, 0.10, 0.30, 0.06, 0.03, 0.48, 0.01],
        T[0.468, 0.468, 0.406, 0.465, 0.464, 0.404, 0.439, 0.476, 0.434, 0.476, 0.421, 0.339, 0.468, 0.006],
        T[0.572, 0.502, 0.785, 0.650, 0.556, 0.780, 0.701, 0.637, 0.640, 0.637, 0.403, 0.201, 0.572, 4.137],
        T[0.571, 0.503, 0.791, 0.656, 0.557, 0.808, 0.740, 0.657, 0.713, 0.657, 0.529, 0.196, 0.571, 4.127],
        T[200, 300, 500, 1000, 1000],
    )
end

function _require_complete(values, label::AbstractString)
    any(ismissing, values) && throw(ArgumentError("$label contains missing values in selected cells"))
    return values
end

"""Construct compact NeuralCrop soil properties from soil codes and pH values."""
function soil_data_from_values(
    soilcode_values::AbstractVector,
    ph_values::AbstractVector,
    selection::CellSelection;
    lookup::SoilLookup{T} = default_soil_lookup(Float32),
    provenance::NamedTuple = (soilcode = nothing, ph = nothing, lookup_version = lookup.version),
) where {T <: AbstractFloat}
    length(soilcode_values) == length(selection.cell_ids) || throw(DimensionMismatch("soilcode count must match selection"))
    length(ph_values) == length(selection.cell_ids) || throw(DimensionMismatch("pH count must match selection"))
    _require_complete(soilcode_values, "soilcode")
    _require_complete(ph_values, "soil pH")

    soilcode = Int32.(soilcode_values)
    class_count = length(lookup.sand)
    all(code -> 1 <= code <= class_count, soilcode) ||
        throw(ArgumentError("soilcode must be in 1:$class_count"))
    ph = T.(ph_values)
    all(value -> isfinite(value) && 0 < value <= 14, ph) ||
        throw(ArgumentError("soil pH must be finite and in (0, 14]"))

    indices = Int.(soilcode)
    saturation_cell = lookup.saturation[indices]
    saturation = repeat(reshape(saturation_cell, 1, :), length(lookup.layer_depth), 1)
    return SoilData(
        selection,
        soilcode,
        ph,
        saturation,
        lookup.sand[indices],
        lookup.silt[indices],
        lookup.clay[indices],
        lookup.diffusivity_dry[indices],
        lookup.diffusivity_15[indices],
        copy(lookup.layer_depth),
        provenance,
    )
end

"""Read aligned soil-code and pH grids and apply the current soil lookup."""
function read_soil_data(
    soilcode_spec::DatasetSpec,
    ph_spec::DatasetSpec,
    grid::GridIndex;
    selection::CellSelection = all_cells(grid),
    lookup::SoilLookup{T} = default_soil_lookup(Float32),
) where {T <: AbstractFloat}
    soilcode = read_static_cell(soilcode_spec, grid; selection)
    ph = read_static_cell(ph_spec, grid; selection)
    provenance = (
        soilcode = soilcode.provenance,
        ph = ph.provenance,
        lookup_version = lookup.version,
    )
    return soil_data_from_values(soilcode.values, ph.values, selection; lookup, provenance)
end

function read_soil_data(
    catalog::DatasetCatalog,
    grid::GridIndex;
    kwargs...,
)
    return read_soil_data(dataset(catalog, :soilcode), dataset(catalog, :soilph), grid; kwargs...)
end

"""Return the soil-property NamedTuple consumed by NeuralCrop's current loader."""
function soil_properties(soil::SoilData)
    return (
        soilcode = soil.soilcode,
        soilph = soil.ph,
        w_sat = soil.saturation,
        sand = soil.sand,
        silt = soil.silt,
        clay = soil.clay,
        tdiff_0 = soil.diffusivity_dry,
        tdiff_15 = soil.diffusivity_15,
        soildepth = soil.layer_depth,
    )
end
