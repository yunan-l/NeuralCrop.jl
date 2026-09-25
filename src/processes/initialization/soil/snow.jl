"""Snow storage, phase-change fluxes, and surface snow diagnostics."""
mutable struct SoilSnow{A}
    pack::A        # Snow water-equivalent storage (mm).
    melt::A        # Snowmelt released today (mm day⁻¹).
    sublimation::A # Snow sublimation loss today (mm day⁻¹).
    runoff::A      # Snow water bypass/runoff loss today (mm day⁻¹).
    height::A      # Physical snow depth (m).
    fraction::A    # Fraction of ground covered by snow (0–1).
end

init_soil_snow(cell_size::Int, device) = init_soil_snow(Float32, cell_size, device)
function init_soil_snow(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    cell_state() = device(zeros(T, cell_size))
    return SoilSnow(ntuple(_ -> cell_state(), 6)...)
end
