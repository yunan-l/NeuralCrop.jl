function _source_units(variable)
    return haskey(variable.attrib, "units") ? String(variable.attrib["units"]) : ""
end

function _validate_units(spec::DatasetSpec, variable)
    source_units = _source_units(variable)
    if !isempty(spec.units) && source_units != spec.units
        throw(ArgumentError("$(spec.variable) units '$source_units' do not match expected units '$(spec.units)'"))
    end
    return source_units
end

function _validate_coordinates(ds, source_names, canonical_names, grid::GridIndex)
    for (canonical_name, expected) in ((:longitude, grid.longitude), (:latitude, grid.latitude))
        position = _dimension_position(canonical_names, canonical_name)
        coordinate_name = String(source_names[position])
        haskey(ds, coordinate_name) || throw(ArgumentError("missing coordinate variable '$coordinate_name'"))
        coordinate = _read_all(ds[coordinate_name])
        length(coordinate) == length(expected) || throw(DimensionMismatch("$canonical_name coordinate length mismatch"))
        all(isapprox.(coordinate, expected; atol = 1.0e-6, rtol = 0)) ||
            throw(ArgumentError("$canonical_name coordinates are not aligned to the canonical grid"))
    end
    return nothing
end

_preserving_selector(selector::Integer) = selector:selector
_preserving_selector(selector) = selector

function _reorder_compact(values, dimensions::Tuple, order::Tuple)
    Set(dimensions) == Set(order) || throw(ArgumentError("requested dimensions $order do not match available dimensions $dimensions"))
    permutation = ntuple(index -> findfirst(==(order[index]), dimensions), length(order))
    return permutedims(values, permutation)
end

"""Read one NetCDF variable, align it to `grid`, and compact space to `:cell`."""
function read_compact_variable(
    spec::DatasetSpec,
    grid::GridIndex;
    selection::CellSelection = all_cells(grid),
    selectors::NamedTuple = NamedTuple(),
    order::Union{Nothing, Tuple} = nothing,
    T::Union{Nothing, Type{<:Number}} = nothing,
)
    return NCDataset(spec.path, "r") do ds
        haskey(ds, spec.variable) || throw(ArgumentError("variable '$(spec.variable)' not found in $(spec.path)"))
        variable = ds[spec.variable]
        source_names = Symbol.(dimnames(variable))
        canonical_names = _canonical_dimension.(source_names)
        allunique(canonical_names) || throw(ArgumentError("variable dimensions are ambiguous after canonicalization"))
        all(name -> name in canonical_names, keys(selectors)) ||
            throw(ArgumentError("a requested selector is not present in $(spec.variable)"))
        longitude_position = _dimension_position(canonical_names, :longitude)
        latitude_position = _dimension_position(canonical_names, :latitude)
        _validate_coordinates(ds, source_names, canonical_names, grid)
        source_units = _validate_units(spec, variable)

        indices = ntuple(ndims(variable)) do position
            name = canonical_names[position]
            selector = haskey(selectors, name) ? getproperty(selectors, name) : Colon()
            return _preserving_selector(selector)
        end
        raw = variable[indices...]
        values = compact_spatial(raw, grid, longitude_position, latitude_position; selection)
        nonspatial_names = Tuple(canonical_names[position] for position in eachindex(canonical_names)
                                 if position ∉ (longitude_position, latitude_position))
        dimensions = (nonspatial_names..., :cell)
        if !isnothing(order)
            values = _reorder_compact(values, dimensions, order)
            dimensions = order
        end
        if !isnothing(T)
            values = T.(values)
        end
        provenance = DataProvenance(DATA_SCHEMA_VERSION, abspath(spec.path), spec.variable, source_units)
        return CompactVariable(values, dimensions, selection, provenance)
    end
end

"""Read a static longitude/latitude variable as one canonical cell vector."""
function read_static_cell(
    spec::DatasetSpec,
    grid::GridIndex;
    selection::CellSelection = all_cells(grid),
    T::Union{Nothing, Type{<:Number}} = nothing,
)
    compact = read_compact_variable(spec, grid; selection, order = (:cell,), T)
    return compact
end
