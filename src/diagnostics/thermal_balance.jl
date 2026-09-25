"""Daily diagnostics for soil enthalpy, ice, and freeze-thaw state."""
mutable struct ThermalBalance{M <: AbstractArray{<:AbstractFloat}}
    surface_energy_flux::M            # Net energy entering soil surface (J m⁻² day⁻¹).
    energy_residual::M                # Daily soil-column energy closure error (J m⁻²).
    untracked_water_energy_flux::M    # Energy correction for externally changed water (J m⁻² day⁻¹).
    rain_energy_input::M              # Sensible heat delivered by rainfall (J m⁻² day⁻¹).
    snowmelt_energy_input::M          # Enthalpy delivered by snowmelt (J m⁻² day⁻¹).
    lateral_runoff_energy_output::M   # Enthalpy removed by lateral runoff (J m⁻² day⁻¹).
    bottom_drainage_energy_output::M  # Enthalpy removed by bottom drainage (J m⁻² day⁻¹).
    percolation_energy_residual::M    # Numerical closure error of flow-energy routing (J m⁻²).
    column_energy::M                  # Layer-integrated soil-column enthalpy (J m⁻²).
    total_ice_storage::M              # Total soil ice (mm water equivalent).
    wilting_ice_storage::M            # Ice within wilting-point water (mm water equivalent).
    available_ice_storage::M          # Ice within plant-available water (mm water equivalent).
    free_ice_storage::M               # Ice within gravitational water (mm water equivalent).
    ice_pool_residual::M              # Difference between total and component ice pools (mm).
    maximum_frozen_fraction::M        # Maximum frozen-water fraction among soil layers (0–1).
    minimum_temperature::M            # Minimum soil-layer temperature (°C).
    maximum_temperature::M            # Maximum soil-layer temperature (°C).
end

function init_thermal_balance(number_of_days::Integer,
                              number_of_cells::Integer,
                              device = identity;
                              T::Type{<:AbstractFloat} = Float32)
    allocate() = device(zeros(T, number_of_days, number_of_cells))
    return ThermalBalance(
        allocate(), allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(), allocate(), allocate(),
    )
end

function record_thermal_balance!(thermal_balance::ThermalBalance,
                                 day_index::Integer,
                                 soil)
    @views begin
        thermal_balance.surface_energy_flux[day_index, :] .=
            soil_thermal_fluxes(soil).surface_energy_flux
        thermal_balance.energy_residual[day_index, :] .=
            soil_thermal_fluxes(soil).energy_residual
        thermal_balance.untracked_water_energy_flux[day_index, :] .=
            soil_thermal_fluxes(soil).untracked_water_energy_flux
        thermal_balance.rain_energy_input[day_index, :] .=
            soil_thermal_fluxes(soil).rain_energy_input
        thermal_balance.snowmelt_energy_input[day_index, :] .=
            soil_thermal_fluxes(soil).snowmelt_energy_input
        thermal_balance.lateral_runoff_energy_output[day_index, :] .=
            soil_thermal_fluxes(soil).lateral_runoff_energy_output
        thermal_balance.bottom_drainage_energy_output[day_index, :] .=
            soil_thermal_fluxes(soil).bottom_drainage_energy_output
        thermal_balance.percolation_energy_residual[day_index, :] .=
            soil_thermal_fluxes(soil).percolation_energy_residual
        thermal_balance.column_energy[day_index, :] .= vec(sum(
            soil_thermal_prognostic(soil).enthalpy .* reshape(
                soil_properties(soil).layer_depth .* eltype(soil_properties(soil).layer_depth)(0.001),
                :, 1,
            ); dims = 1,
        ))
        thermal_balance.total_ice_storage[day_index, :] .=
            vec(sum(soil_water_prognostic(soil).ice_storage; dims = 1))
        thermal_balance.wilting_ice_storage[day_index, :] .= vec(sum(
            soil_water_prognostic(soil).wilting_ice_fraction .* soil_water_auxiliary(soil).wilting_storage;
            dims = 1,
        ))
        thermal_balance.available_ice_storage[day_index, :] .=
            vec(sum(soil_water_prognostic(soil).available_ice_storage; dims = 1))
        thermal_balance.free_ice_storage[day_index, :] .=
            vec(sum(soil_water_prognostic(soil).free_ice_storage; dims = 1))
        component_ice =
            soil_water_prognostic(soil).wilting_ice_fraction .* soil_water_auxiliary(soil).wilting_storage .+
            soil_water_prognostic(soil).available_ice_storage .+
            soil_water_prognostic(soil).free_ice_storage
        thermal_balance.ice_pool_residual[day_index, :] .=
            vec(maximum(abs.(soil_water_prognostic(soil).ice_storage .- component_ice); dims = 1))
        thermal_balance.maximum_frozen_fraction[day_index, :] .=
            vec(maximum(soil_thermal_prognostic(soil).frozen_fraction; dims = 1))
        thermal_balance.minimum_temperature[day_index, :] .=
            vec(minimum(soil_thermal_prognostic(soil).temperature; dims = 1))
        thermal_balance.maximum_temperature[day_index, :] .=
            vec(maximum(soil_thermal_prognostic(soil).temperature; dims = 1))
    end
    return nothing
end
