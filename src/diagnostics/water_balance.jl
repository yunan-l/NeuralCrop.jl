"""
    WaterBalance

Optional daily water-budget diagnostics. All fields have shape
`(number_of_days, number_of_cells)` and use millimetres of water.

`residual` is positive when water enters the model but is neither retained in
soil, snow, or surface-litter storage nor represented by a recorded outgoing
flux. Surface-litter interception is an internal transfer; litter evaporation
is an outgoing flux.
"""
mutable struct WaterBalance{M <: AbstractArray{<:AbstractFloat}}
    precipitation::M            # Atmospheric water input (mm day⁻¹).
    rain_after_snow::M          # Liquid precipitation remaining after snow processing (mm day⁻¹).
    soil_storage_before::M      # Liquid soil-water stock at start of day (mm).
    soil_storage_after::M       # Liquid soil-water stock at end of day (mm).
    soil_ice_storage_before::M  # Soil-ice stock at start of day (mm water equivalent).
    soil_ice_storage_after::M   # Soil-ice stock at end of day (mm water equivalent).
    snow_storage_before::M      # Snow stock at start of day (mm water equivalent).
    snow_storage_after::M       # Snow stock at end of day (mm water equivalent).
    litter_storage_before::M    # Surface-litter water stock at start of day (mm).
    litter_storage_after::M     # Surface-litter water stock at end of day (mm).
    snowmelt::M                 # Water released by snow melt (mm day⁻¹).
    snow_sublimation::M         # Snow sublimation loss (mm day⁻¹).
    snow_runoff::M              # Snow bypass/runoff loss (mm day⁻¹).
    unaccounted_snow_flux::M    # Snow-subsystem water closure error (mm).
    interception::M             # Rain intercepted and evaporated by crop canopy (mm day⁻¹).
    litter_interception::M      # Rain transferred into surface-litter storage (mm day⁻¹).
    litter_evaporation::M       # Evaporation from surface litter (mm day⁻¹).
    transpiration::M            # Crop transpiration loss (mm day⁻¹).
    evaporation::M              # Bare-soil evaporation loss (mm day⁻¹).
    surface_runoff::M           # Surface runoff loss (mm day⁻¹).
    lateral_runoff::M           # Layer-integrated lateral runoff loss (mm day⁻¹).
    bottom_drainage::M          # Drainage loss from bottom layer (mm day⁻¹).
    remaining_infiltration::M   # Unprocessed infiltration retained in the ledger (mm day⁻¹).
    residual::M                 # Absolute daily water-budget closure error (mm).
end

"""
    init_water_balance(number_of_days, number_of_cells, device=identity; T=Float32)

Allocate an optional daily water-balance ledger on the selected device.
Currently the ledger supports rainfed simulations only.
"""
function init_water_balance(number_of_days::Integer,
                            number_of_cells::Integer,
                            device = identity;
                            T::Type{<:AbstractFloat} = Float32)
    allocate() = device(zeros(T, number_of_days, number_of_cells))

    return WaterBalance(
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(),
    )
end

function record_water_balance_start!(water_balance::WaterBalance,
                                     day_index::Integer,
                                     soil,
                                     precipitation)
    @views water_balance.precipitation[day_index, :] .= precipitation
    @views water_balance.soil_storage_before[day_index, :] .= vec(sum(soil_water_prognostic(soil).storage; dims = 1))
    @views water_balance.soil_ice_storage_before[day_index, :] .=
        vec(sum(soil_water_prognostic(soil).ice_storage; dims = 1))
    @views water_balance.snow_storage_before[day_index, :] .= soil_snow_prognostic(soil).pack
    @views water_balance.litter_storage_before[day_index, :] .=
        soil_surface_litter_prognostic(soil).water_storage
    return nothing
end

function record_water_balance_after_snow!(water_balance::WaterBalance,
                                          day_index::Integer,
                                          precipitation)
    @views water_balance.rain_after_snow[day_index, :] .= precipitation
    return nothing
end

function record_water_balance_end!(water_balance::WaterBalance,
                                   day_index::Integer,
                                   soil,
                                   crop)
    @views begin
        water_balance.soil_storage_after[day_index, :] .= vec(sum(soil_water_prognostic(soil).storage; dims = 1))
        water_balance.soil_ice_storage_after[day_index, :] .=
            vec(sum(soil_water_prognostic(soil).ice_storage; dims = 1))
        water_balance.snow_storage_after[day_index, :] .= soil_snow_prognostic(soil).pack
        water_balance.litter_storage_after[day_index, :] .=
            soil_surface_litter_prognostic(soil).water_storage
        water_balance.snowmelt[day_index, :] .= soil_snow_fluxes(soil).melt
        water_balance.snow_sublimation[day_index, :] .= soil_snow_fluxes(soil).sublimation
        water_balance.snow_runoff[day_index, :] .= soil_snow_fluxes(soil).runoff
        water_balance.interception[day_index, :] .= crop_fluxes(crop).water.interception
        water_balance.litter_interception[day_index, :] .=
            soil_surface_litter_fluxes(soil).interception
        water_balance.litter_evaporation[day_index, :] .=
            soil_surface_litter_fluxes(soil).evaporation
        water_balance.transpiration[day_index, :] .= vec(sum(crop_fluxes(crop).water.transpiration_layer; dims = 1))
        water_balance.evaporation[day_index, :] .= vec(sum(soil_water_fluxes(soil).evaporation; dims = 1))
        water_balance.surface_runoff[day_index, :] .= soil_water_fluxes(soil).surface_runoff
        water_balance.lateral_runoff[day_index, :] .= vec(sum(soil_water_fluxes(soil).lateral_runoff; dims = 1))
        water_balance.bottom_drainage[day_index, :] .= soil_water_fluxes(soil).bottom_drainage
        water_balance.remaining_infiltration[day_index, :] .= soil_water_fluxes(soil).infiltration

        water_balance.unaccounted_snow_flux[day_index, :] .=
            water_balance.snow_storage_before[day_index, :] .+
            water_balance.precipitation[day_index, :] .-
            water_balance.rain_after_snow[day_index, :] .-
            water_balance.snow_storage_after[day_index, :] .-
            water_balance.snow_sublimation[day_index, :] .-
            water_balance.snow_runoff[day_index, :]

        water_balance.residual[day_index, :] .=
            water_balance.soil_storage_before[day_index, :] .+
            water_balance.soil_ice_storage_before[day_index, :] .+
            water_balance.snow_storage_before[day_index, :] .+
            water_balance.litter_storage_before[day_index, :] .+
            water_balance.precipitation[day_index, :] .-
            water_balance.soil_storage_after[day_index, :] .-
            water_balance.soil_ice_storage_after[day_index, :] .-
            water_balance.snow_storage_after[day_index, :] .-
            water_balance.litter_storage_after[day_index, :] .-
            water_balance.interception[day_index, :] .-
            water_balance.litter_evaporation[day_index, :] .-
            water_balance.transpiration[day_index, :] .-
            water_balance.evaporation[day_index, :] .-
            water_balance.snow_sublimation[day_index, :] .-
            water_balance.snow_runoff[day_index, :] .-
            water_balance.surface_runoff[day_index, :] .-
            water_balance.lateral_runoff[day_index, :] .-
            water_balance.bottom_drainage[day_index, :]
    end

    return nothing
end
