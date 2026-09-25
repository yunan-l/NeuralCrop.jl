"""Hydrological and thermal state of the above-ground litter layer."""
mutable struct SoilSurfaceLitter{A}
    dry_matter::A    # Above-ground litter dry-matter stock (g dry matter m⁻²).
    depth::A         # Effective surface-litter layer thickness (m).
    cover::A         # Fractional ground cover by above-ground litter (0–1).
    water_capacity::A # Maximum litter water-holding capacity (mm).
    water_storage::A # Liquid water stored in surface litter (mm).
    interception::A  # Rainfall intercepted by litter today (mm day⁻¹).
    evaporation::A   # Water evaporated from surface litter today (mm day⁻¹).
    temperature::A   # Surface-litter temperature (°C).
    conductivity::A  # Effective litter thermal conductivity (W m⁻¹ K⁻¹).
end

init_soil_surface_litter(cell_size::Int, device) =
    init_soil_surface_litter(Float32, cell_size, device)
function init_soil_surface_litter(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    state() = device(zeros(T, cell_size))
    return SoilSurfaceLitter(ntuple(_ -> state(), 9)...)
end
