"""Static or slowly varying soil physical and chemical properties."""
mutable struct SoilProperties{A, M}
    sand_fraction::M # Soil sand mass fraction by layer (0–1).
    clay_fraction::M # Soil clay mass fraction by layer (0–1).
    ph::A            # Soil pH used by nitrogen transformations.
    layer_depth::A   # Thickness of each model soil layer (mm).
end

init_soil_properties(cell_size::Int, soildepth, device) =
    init_soil_properties(Float32, cell_size, soildepth, device)
function init_soil_properties(::Type{T}, cell_size::Int, soildepth, device) where {T <: AbstractFloat}
    return SoilProperties(
        device(zeros(T, 1, cell_size)),
        device(zeros(T, 1, cell_size)),
        device(zeros(T, cell_size)),
        device(T.(soildepth)),
    )
end
