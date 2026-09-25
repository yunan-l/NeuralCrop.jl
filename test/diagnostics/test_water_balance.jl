using NeuralCrop
using Test

@testset "Daily water-balance diagnostics" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)
    state = test_model_state(crop, soil)
    water_balance = init_water_balance(3, 1, identity)

    soil.water.storage .= 20.0f0
    NeuralCrop.record_water_balance_start!(water_balance, 1, state, Float32[10.0])
    NeuralCrop.record_water_balance_after_snow!(water_balance, 1, Float32[10.0])

    soil.water.storage[1, 1] = 24.5f0
    crop.fluxes.water.interception .= 1.0f0
    crop.fluxes.water.transpiration_layer .= 0.4f0
    soil.water.evaporation .= 0.2f0
    soil.water.surface_runoff .= 0.5f0
    soil.water.lateral_runoff .= 0.1f0
    soil.water.bottom_drainage .= 0.5f0
    NeuralCrop.record_water_balance_end!(water_balance, 1, state, state)

    @test water_balance.soil_storage_before[1, 1] == 100.0f0
    @test water_balance.soil_storage_after[1, 1] == 104.5f0
    @test water_balance.transpiration[1, 1] == 2.0f0
    @test water_balance.evaporation[1, 1] == 1.0f0
    @test water_balance.lateral_runoff[1, 1] == 0.5f0
    @test water_balance.unaccounted_snow_flux[1, 1] == 0.0f0
    @test water_balance.residual[1, 1] ≈ 0.0f0 atol = 1.0f-6

    soil.water.storage .= 20.0f0
    soil.snow.pack .= 2.0f0
    NeuralCrop.record_water_balance_start!(water_balance, 2, state, Float32[3.0])
    NeuralCrop.record_water_balance_after_snow!(water_balance, 2, Float32[0.0])
    soil.snow.pack .= 4.9f0
    soil.snow.sublimation .= 0.1f0
    crop.fluxes.water.interception .= 0.0f0
    crop.fluxes.water.transpiration_layer .= 0.0f0
    soil.water.evaporation .= 0.0f0
    soil.water.surface_runoff .= 0.0f0
    soil.water.lateral_runoff .= 0.0f0
    soil.water.bottom_drainage .= 0.0f0
    NeuralCrop.record_water_balance_end!(water_balance, 2, state, state)

    @test water_balance.snow_sublimation[2, 1] == 0.1f0
    @test water_balance.unaccounted_snow_flux[2, 1] ≈ 0.0f0 atol = 1.0f-6
    @test water_balance.residual[2, 1] ≈ 0.0f0 atol = 1.0f-6

    # Litter interception is an internal transfer. Its storage change and
    # subsequent evaporation must nevertheless close the total water budget.
    soil.snow.pack .= 0.0f0
    soil.snow.sublimation .= 0.0f0
    soil.water.storage .= 20.0f0
    soil.surface_litter.water_storage .= 1.0f0
    NeuralCrop.record_water_balance_start!(water_balance, 3, state, Float32[4.0])
    NeuralCrop.record_water_balance_after_snow!(water_balance, 3, Float32[4.0])
    soil.surface_litter.interception .= 2.0f0
    soil.surface_litter.evaporation .= 0.5f0
    soil.surface_litter.water_storage .= 2.5f0
    soil.water.storage[1, 1] = 22.0f0
    crop.fluxes.water.interception .= 0.0f0
    crop.fluxes.water.transpiration_layer .= 0.0f0
    soil.water.evaporation .= 0.0f0
    soil.water.surface_runoff .= 0.0f0
    soil.water.lateral_runoff .= 0.0f0
    soil.water.bottom_drainage .= 0.0f0
    NeuralCrop.record_water_balance_end!(water_balance, 3, state, state)

    @test water_balance.litter_interception[3, 1] == 2.0f0
    @test water_balance.litter_evaporation[3, 1] == 0.5f0
    @test water_balance.residual[3, 1] ≈ 0.0f0 atol = 1.0f-6
end
