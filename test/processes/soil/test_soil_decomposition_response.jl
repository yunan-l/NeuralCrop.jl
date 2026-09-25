using NeuralCrop
using Test

@testset "LPJmL soil decomposition response and rates" begin
    @test lpjmlparams.k_soil10.fast ≈ 0.04f0 / 365.0f0
    @test lpjmlparams.k_soil10.slow ≈ 0.001f0 / 365.0f0
    @test lpjmlparams.atmfrac == 0.5f0

    soil = init_soil(2, soilparams.soildepth, identity)
    state = test_model_state(soil)
    soil.thermal.temperature .= 10.0f0
    soil.water.saturation_storage .= 100.0f0
    soil.water.holding_capacity_storage .= 100.0f0
    soil.water.wilting_storage .= 0.0f0
    soil.water.free_water .= 0.0f0
    soil.water.relative_content[:, 1] .= 0.5f0
    soil.water.relative_content[:, 2] .= 0.25f0
    soil.water.available_ice_storage[:, 2] .= 50.0f0
    soil.surface_litter.temperature .= -20.0f0
    soil.surface_litter.water_capacity .= 1.0f0
    soil.surface_litter.water_storage .= 0.5f0

    soil_decomp_response!(state)

    # Cell 2 has half of its pore volume occupied by ice. Its 25 mm liquid
    # therefore has the same 50% liquid-pore saturation as cell 1's 50 mm.
    @test soil.decomposition.response[:, 1] ≈
          soil.decomposition.response[:, 2] rtol = 2.0f-6
    @test all(iszero, soil.decomposition.litter_response[1, :])
    @test soil.decomposition.litter_response[2, :] ≈
          soil.decomposition.response[1, :]
    @test soil.decomposition.litter_response[3, :] ≈
          soil.decomposition.response[1, :]

    crop = init_crop(1, identity)
    decay_soil = init_soil(1, soilparams.soildepth, identity)
    decay_state = test_model_state(crop, decay_soil)
    decay_soil.thermal.temperature .= 10.0f0
    decay_soil.water.saturation_storage .= 100.0f0
    decay_soil.water.holding_capacity_storage .= 100.0f0
    decay_soil.water.relative_content .= 0.5f0
    decay_soil.carbon.fast .= 100.0f0
    decay_soil.carbon.slow .= 100.0f0
    soil_carbon!(decay_state, decay_state)
    response = decay_soil.decomposition.response[1, 1]
    @test decay_soil.carbon.decomposed_fast[1, 1] ≈
        -100.0f0 * expm1(-0.04f0 / 365.0f0 * response) rtol = 2.0f-5
    @test decay_soil.carbon.decomposed_slow[1, 1] ≈
        -100.0f0 * expm1(-0.001f0 / 365.0f0 * response) rtol = 2.0f-4
end
