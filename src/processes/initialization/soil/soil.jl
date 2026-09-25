"""Process-grouped soil state with CPU/GPU-compatible array leaves."""
mutable struct Soil{P, W, T, C, N, D, G, L, S}
    properties::P     # Static soil texture, pH, and layer geometry.
    water::W          # Layer water/ice stocks, hydraulic properties, and daily water fluxes.
    thermal::T        # Layer temperature, enthalpy, phase state, and energy fluxes.
    carbon::C         # Litter and soil organic-carbon pools and decomposition fluxes.
    nitrogen::N       # Mineral/organic nitrogen pools and transformation fluxes.
    decomposition::D  # Environmental responses, fixed C/N routing configuration, and workspace.
    management::G     # Tillage and bioturbation routing/diagnostics.
    surface_litter::L # Above-ground litter hydrological and thermal state.
    snow::S           # Snow water stock and current-day snow fluxes.
end

init_soil(cell_size::Int, soildepth::AbstractArray{T}, device; kwargs...) where {T <: AbstractFloat} =
    init_soil(T, cell_size, soildepth, device; kwargs...)
function init_soil(::Type{T},
                   cell_size::Int,
                   soildepth::AbstractArray{S},
                   device;
                   litc_layers::Int = 3,
                   soil_layers::Int = 5) where {T <: AbstractFloat, S <: AbstractFloat}
    return Soil(
        init_soil_properties(T, cell_size, soildepth, device),
        init_soil_water(T, cell_size, device; soil_layers = soil_layers),
        init_soil_thermal(T, cell_size, device; soil_layers = soil_layers),
        init_soil_carbon(T, cell_size, device;
                         litter_layers = litc_layers, soil_layers = soil_layers),
        init_soil_nitrogen(T, cell_size, device;
                           litter_layers = litc_layers, soil_layers = soil_layers),
        init_soil_decomposition(T, cell_size, device; soil_layers = soil_layers),
        init_soil_management(T, cell_size, device; litter_layers = litc_layers),
        init_soil_surface_litter(T, cell_size, device),
        init_soil_snow(T, cell_size, device),
    )
end
