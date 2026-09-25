"""Crop-management inputs and grid-cell location state."""
mutable struct ManagedLand{A}
    manure::A           # Prescribed manure application (gN m⁻² crop⁻¹).
    fertilizer::A       # Prescribed mineral fertilizer application (gN m⁻² crop⁻¹).
    residue_fraction::A # Fraction of harvested above-ground residue retained (0–1).
    latitude::A         # Grid-cell latitude (degrees north).
end

init_managed_land(cell_size::Int, device) = init_managed_land(Float32, cell_size, device)
function init_managed_land(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    cell_state() = device(zeros(T, cell_size))
    return ManagedLand(ntuple(_ -> cell_state(), 4)...)
end
