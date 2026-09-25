const SURFACE_LITTER = 1
const INCORPORATED_LITTER = 2
const ROOT_LITTER = 3

"""Soil management operators and daily internal litter-routing diagnostics."""
mutable struct SoilManagement{M, A}
    tillage_fraction::M    # Routing matrix among litter classes during cultivation (fraction).
    tillage_density_factor::M # Tilled-layer bulk density relative to settled soil (0–1; prognostic).
    tillage_carbon::A      # Carbon redistributed by tillage today (gC m⁻² day⁻¹).
    tillage_nitrogen::A    # Nitrogen redistributed by tillage today (gN m⁻² day⁻¹).
    bioturbation_carbon::A # Carbon redistributed vertically by bioturbation today (gC m⁻² day⁻¹).
    bioturbation_nitrogen::A # Nitrogen redistributed vertically by bioturbation today (gN m⁻² day⁻¹).
end

init_soil_management(cell_size::Int, device; kwargs...) =
    init_soil_management(Float32, cell_size, device; kwargs...)
function init_soil_management(::Type{T}, cell_size::Int, device;
                              litter_layers::Int = 3,
                              tillage_layers::Int = 1) where {T <: AbstractFloat}
    cell_state() = device(zeros(T, cell_size))
    return SoilManagement(
        device(zeros(T, litter_layers, litter_layers)),
        device(ones(T, tillage_layers, cell_size)),
        cell_state(), cell_state(), cell_state(), cell_state(),
    )
end
