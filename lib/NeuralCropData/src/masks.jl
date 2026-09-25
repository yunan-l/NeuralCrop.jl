function _valid_crop_fraction(value, threshold)
    ismissing(value) && return false
    value isa AbstractFloat && !isfinite(value) && throw(ArgumentError("landuse contains a non-finite value"))
    return value > threshold
end

"""
Build a fixed allocation selection from the temporal union of positive land use.
The returned `active` matrix retains annual activity only for allocated cells.
"""
function build_crop_mask(
    grid::GridIndex,
    landuse::AbstractMatrix;
    selection::CellSelection = all_cells(grid),
    threshold::Real = 0,
)
    size(landuse, 2) == length(selection.cell_ids) ||
        throw(DimensionMismatch("landuse must have dimensions time × selected cell"))
    active_full = BitMatrix(undef, size(landuse))
    for index in eachindex(landuse)
        active_full[index] = _valid_crop_fraction(landuse[index], threshold)
    end
    allocated_local = findall(vec(any(active_full; dims = 1)))
    isempty(allocated_local) && throw(ArgumentError("selected CFT has no active land-use cells"))
    compact_indices = selection.compact_indices[allocated_local]
    allocated = select_cells(grid, compact_indices)

    value_type = Base.nonmissingtype(eltype(landuse))
    fraction = Matrix{value_type}(undef, size(landuse, 1), length(allocated_local))
    for (output_cell, input_cell) in pairs(allocated_local), time_index in axes(landuse, 1)
        value = landuse[time_index, input_cell]
        fraction[time_index, output_cell] = ismissing(value) ? zero(value_type) : value_type(value)
    end
    return CropMask(allocated, fraction, active_full[:, allocated_local])
end

"""Build patches for one CFT and water regime from its positive land-use cells."""
function build_patch_domain(
    grid::GridIndex,
    landuse::AbstractMatrix,
    cft_id::Integer;
    irrigated::Bool = false,
    selection::CellSelection = all_cells(grid),
    valid::Union{Nothing, AbstractVector{Bool}} = nothing,
    threshold::Real = 0,
)
    mask = build_crop_mask(grid, landuse; selection, threshold)
    valid_cells = isnothing(valid) ? trues(length(mask.selection.cell_ids)) : BitVector(valid)
    length(valid_cells) == length(mask.selection.cell_ids) || throw(DimensionMismatch(
        "valid must contain one value per land-use-selected cell",
    ))
    any(valid_cells) || throw(ArgumentError("no valid patch cells remain"))
    allocated = mask.selection
    patch_selection = CellSelection(
        allocated.compact_indices[valid_cells], allocated.cell_ids[valid_cells],
    )
    landfrac = vec(maximum(mask.fraction[:, valid_cells]; dims = 1))
    patch_count = length(patch_selection.cell_ids)
    return PatchDomain(
        1:patch_count,
        patch_selection.compact_indices,
        patch_selection.cell_ids,
        fill(Int32(cft_id), patch_count),
        fill(irrigated, patch_count),
        landfrac,
    )
end

"""Combine CFT/water-regime patch domains while preserving duplicate grid cells."""
function combine_patch_domains(domains::AbstractVector{<:PatchDomain})
    isempty(domains) && throw(ArgumentError("at least one patch domain is required"))
    T = promote_type((eltype(domain.landfrac) for domain in domains)...)
    compact_indices = reduce(vcat, (domain.compact_indices for domain in domains))
    cell_ids = reduce(vcat, (domain.cell_ids for domain in domains))
    cft_ids = reduce(vcat, (domain.cft_ids for domain in domains))
    irrigated = reduce(vcat, (domain.irrigated for domain in domains))
    landfrac = reduce(vcat, (T.(domain.landfrac) for domain in domains))
    return PatchDomain(1:length(cell_ids), compact_indices, cell_ids, cft_ids, irrigated, landfrac)
end
