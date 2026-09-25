"""Layered soil temperature, enthalpy, phase state, and thermal properties."""
mutable struct SoilThermal{A, B, M}
    temperature::M                   # Soil-layer temperature (°C).
    enthalpy::M                      # Volumetric soil enthalpy relative to 0 °C (J m⁻³).
    frozen_fraction::M               # Fraction of layer water that is frozen (0–1).
    freeze_depth::M                  # Effective frozen depth within each soil layer (mm).
    heat_capacity_frozen::M          # Frozen-layer volumetric heat capacity (J m⁻³ K⁻¹).
    heat_capacity_unfrozen::M        # Unfrozen-layer volumetric heat capacity (J m⁻³ K⁻¹).
    latent_heat::M                   # Latent heat stored/released by layer phase change (J m⁻³).
    conductivity_frozen::M           # Frozen-layer thermal conductivity (W m⁻¹ K⁻¹).
    conductivity_unfrozen::M         # Unfrozen-layer thermal conductivity (W m⁻¹ K⁻¹).
    water_reference::M               # Water stock used as phase-change reference (mm).
    percolation_energy::M            # Enthalpy carried by layer water flow (J m⁻² day⁻¹).
    surface_energy_flux::A           # Net daily energy entering the soil surface (J m⁻² day⁻¹).
    energy_residual::A               # Soil-column daily energy-balance residual (J m⁻²).
    untracked_water_energy_flux::A   # Energy correction for externally changed water (J m⁻² day⁻¹).
    rain_energy_input::A             # Sensible heat delivered by rainfall (J m⁻² day⁻¹).
    snowmelt_energy_input::A         # Enthalpy delivered by snowmelt (J m⁻² day⁻¹).
    lateral_runoff_energy_output::A  # Enthalpy removed by lateral runoff (J m⁻² day⁻¹).
    bottom_drainage_energy_output::A # Enthalpy removed by bottom drainage (J m⁻² day⁻¹).
    percolation_energy_residual::A   # Numerical closure residual of flow-energy routing (J m⁻²).
    initialized::B                   # Thermal-profile initialization flag per grid cell.
    diffusivity_0::A                 # Soil thermal diffusivity at zero water content (mm² s⁻¹).
    diffusivity_15::A                # Soil thermal diffusivity at reference moisture 0.15 (mm² s⁻¹).
end

init_soil_thermal(cell_size::Int, device; kwargs...) =
    init_soil_thermal(Float32, cell_size, device; kwargs...)
function init_soil_thermal(::Type{T}, cell_size::Int, device;
                           soil_layers::Int = 5) where {T <: AbstractFloat}
    layer_state() = device(zeros(T, soil_layers, cell_size))
    cell_state() = device(zeros(T, cell_size))
    return SoilThermal(
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        layer_state(),
        cell_state(),
        cell_state(),
        cell_state(),
        cell_state(),
        cell_state(),
        cell_state(),
        cell_state(),
        cell_state(),
        device(zeros(Bool, cell_size)),
        cell_state(),
        cell_state(),
    )
end
