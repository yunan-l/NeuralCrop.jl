using NeuralCrop
using Test

@testset "LPJmL untracked water-mass enthalpy" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)
    state = test_model_state(crop, soil)
    soil.properties.sand_fraction .= 0.4f0
    soil.properties.clay_fraction .= 0.2f0
    soil.water.storage .= reshape(Float32[40, 60, 100, 200, 200], 5, 1)
    pedotransfer!(state)
    soil_temperature!(state, Float32[5.0], Float32[5.0])

    crop.fluxes.water.transpiration_layer .= 0.5f0
    soil.water.evaporation .= 0.25f0
    removed_water = sum(crop.fluxes.water.transpiration_layer) + sum(soil.water.evaporation)
    soil_evapotranspiration!(state, state)
    soil_temperature!(state, Float32[5.0])

    liquid_water_enthalpy = soil_thermal_params.volumetric_fusion_heat +
                            soil_thermal_params.water_heat_capacity * 5.0f0
    expected_energy_change = -removed_water * 0.001f0 * liquid_water_enthalpy
    @test soil.thermal.untracked_water_energy_flux[1] ≈ expected_energy_change rtol = 2.0f-5
    @test maximum(abs.(soil.thermal.temperature .- 5.0f0)) < 2.0f-5
    @test abs(soil.thermal.energy_residual[1]) < 256.0f0
    @test all(isfinite, soil.thermal.enthalpy)
end
