const _LONGITUDE_NAMES = (:longitude, :lon, :x)
const _LATITUDE_NAMES = (:latitude, :lat, :y)
const _TIME_NAMES = (:time, :year)
# Source files may still use LPJmL's historical `pft` dimension name. It is
# normalized at this I/O boundary; all NeuralCropData contracts use `:cft`.
const _CFT_NAMES = (:cft, :pft, :crop, :band)

function _canonical_dimension(name::Union{Symbol, AbstractString})
    symbol = Symbol(lowercase(String(name)))
    symbol in _LONGITUDE_NAMES && return :longitude
    symbol in _LATITUDE_NAMES && return :latitude
    symbol in _TIME_NAMES && return :time
    symbol in _CFT_NAMES && return :cft
    return symbol
end

function _dimension_position(names, target::Symbol)
    positions = findall(==(target), names)
    length(positions) == 1 || throw(ArgumentError("expected exactly one $target dimension, found $(length(positions))"))
    return only(positions)
end

function _read_all(variable)
    indices = ntuple(_ -> Colon(), ndims(variable))
    return variable[indices...]
end

"""Read and validate the canonical `cellid` grid."""
function read_grid(spec::DatasetSpec; T::Type{<:AbstractFloat} = Float64)
    return NCDataset(spec.path, "r") do ds
        haskey(ds, spec.variable) || throw(ArgumentError("variable '$(spec.variable)' not found in $(spec.path)"))
        variable = ds[spec.variable]
        names = _canonical_dimension.(Symbol.(dimnames(variable)))
        ndims(variable) == 2 || throw(ArgumentError("grid cellid variable must be two-dimensional"))
        lon_position = _dimension_position(names, :longitude)
        lat_position = _dimension_position(names, :latitude)

        source_names = Symbol.(dimnames(variable))
        longitude = T.(_read_all(ds[String(source_names[lon_position])]))
        latitude = T.(_read_all(ds[String(source_names[lat_position])]))
        raw = _read_all(variable)
        cellid_source = permutedims(raw, (lon_position, lat_position))

        cellid = fill(Int32(-1), length(longitude), length(latitude))
        records = Tuple{Int32, Int32, Int32}[]
        for latitude_index in eachindex(latitude), longitude_index in eachindex(longitude)
            value = cellid_source[longitude_index, latitude_index]
            if !ismissing(value)
                id = Int32(value)
                id >= 0 || throw(ArgumentError("valid cell ids must be non-negative, found $id"))
                cellid[longitude_index, latitude_index] = id
                push!(records, (id, Int32(longitude_index), Int32(latitude_index)))
            end
        end
        isempty(records) && throw(ArgumentError("grid contains no valid cells"))
        sort!(records; by = first)
        ids = first.(records)
        allunique(ids) || throw(ArgumentError("grid cell ids must be unique"))
        return GridIndex(
            longitude,
            latitude,
            cellid,
            ids,
            getindex.(records, 2),
            getindex.(records, 3),
        )
    end
end

read_grid(path::AbstractString; variable::AbstractString = "cellid", kwargs...) =
    read_grid(DatasetSpec(path, variable); kwargs...)

all_cells(grid::GridIndex) = CellSelection(eachindex(grid.cell_ids), grid.cell_ids)

function _validate_selection(grid::GridIndex, selection::CellSelection)
    all(index -> checkbounds(Bool, grid.cell_ids, index), selection.compact_indices) ||
        throw(BoundsError(grid.cell_ids, selection.compact_indices))
    grid.cell_ids[selection.compact_indices] == selection.cell_ids ||
        throw(ArgumentError("cell selection does not belong to the canonical grid"))
    return selection
end

function select_cells(grid::GridIndex, compact_indices::AbstractVector{<:Integer})
    indices = Int.(compact_indices)
    all(index -> checkbounds(Bool, grid.cell_ids, index), indices) ||
        throw(BoundsError(grid.cell_ids, indices))
    issorted(indices) || throw(ArgumentError("compact cell indices must preserve canonical order"))
    allunique(indices) || throw(ArgumentError("compact cell indices must be unique"))
    return CellSelection(indices, grid.cell_ids[indices])
end

function select_cells(grid::GridIndex, mask::AbstractVector{Bool})
    length(mask) == length(grid.cell_ids) || throw(DimensionMismatch("mask length must match canonical cells"))
    return select_cells(grid, findall(mask))
end

"""
Compact longitude/latitude dimensions into the canonical final `cell` dimension.
All non-spatial dimensions retain their relative order.
"""
function compact_spatial(
    values::AbstractArray,
    grid::GridIndex,
    longitude_dimension::Integer,
    latitude_dimension::Integer;
    selection::CellSelection = all_cells(grid),
)
    _validate_selection(grid, selection)
    ndims(values) >= 2 || throw(ArgumentError("values must contain longitude and latitude dimensions"))
    longitude_dimension != latitude_dimension || throw(ArgumentError("spatial dimensions must be distinct"))
    size(values, longitude_dimension) == length(grid.longitude) ||
        throw(DimensionMismatch("longitude size does not match canonical grid"))
    size(values, latitude_dimension) == length(grid.latitude) ||
        throw(DimensionMismatch("latitude size does not match canonical grid"))

    other_dimensions = [dimension for dimension in 1:ndims(values) if dimension ∉ (longitude_dimension, latitude_dimension)]
    permutation = (other_dimensions..., longitude_dimension, latitude_dimension)
    ordered = permutedims(values, permutation)
    other_sizes = size(ordered)[1:length(other_dimensions)]
    rows = isempty(other_sizes) ? 1 : prod(other_sizes)
    spatial = reshape(ordered, rows, length(grid.longitude) * length(grid.latitude))

    compact_indices = selection.compact_indices
    longitude_indices = Int.(grid.longitude_indices[compact_indices])
    latitude_indices = Int.(grid.latitude_indices[compact_indices])
    linear_indices = longitude_indices .+ (latitude_indices .- 1) .* length(grid.longitude)
    selected = spatial[:, linear_indices]
    return reshape(selected, (other_sizes..., length(selection.cell_ids)))
end

"""Expand one compact cell vector to the canonical longitude/latitude grid."""
function expand_to_grid(
    values::AbstractVector{T},
    grid::GridIndex;
    selection::CellSelection = all_cells(grid),
    fill_value::T = T <: AbstractFloat ? T(NaN) : zero(T),
) where {T}
    _validate_selection(grid, selection)
    length(values) == length(selection.cell_ids) || throw(DimensionMismatch("value count must match selection"))
    output = fill(fill_value, length(grid.longitude), length(grid.latitude))
    for (local_index, compact_index) in pairs(selection.compact_indices)
        longitude_index = grid.longitude_indices[compact_index]
        latitude_index = grid.latitude_indices[compact_index]
        output[longitude_index, latitude_index] = values[local_index]
    end
    return output
end
