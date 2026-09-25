using NeuralCrop
using Test

@testset "LPJmL daily soil decomposition flux guards" begin
    crop = init_crop(2, identity)
    soil = init_soil(2, soilparams.soildepth, identity)
    state = test_model_state(crop, soil)

    soil.carbon.litter .= 10.0f0
    soil.nitrogen.litter .= 1.0f0
    soil.carbon.fast .= 10.0f0
    soil.carbon.slow .= 10.0f0
    soil.nitrogen.fast .= 1.0f0
    soil.nitrogen.slow .= 1.0f0
    soil.carbon.litter_response .= 1.0f0
    soil.nitrogen.litter_response .= 1.0f0
    soil.surface_litter.temperature .= 10.0f0
    soil.water.saturation_storage .= 100.0f0
    soil.water.holding_capacity_storage .= 100.0f0
    soil.water.relative_content .= 0.5f0
    soil.thermal.temperature[:, 1] .= -20.0f0
    soil.thermal.temperature[:, 2] .= 10.0f0

    soil_carbon!(state, state)
    soil_nitrogen!(state, state)

    # LPJmL gates the entire litter block with top-layer gtemp_soil > 0.
    @test all(iszero, soil.carbon.decomposed_litter[:, 1])
    @test all(iszero, soil.nitrogen.decomposed_litter[:, 1])
    @test all(iszero, soil.carbon.decomposed_fast[:, 1])
    @test all(iszero, soil.carbon.decomposed_slow[:, 1])
    @test all(iszero, soil.nitrogen.decomposed_fast[:, 1])
    @test all(iszero, soil.nitrogen.decomposed_slow[:, 1])
    @test all(soil.carbon.decomposed_litter[:, 2] .> 0.0f0)
    @test all(soil.nitrogen.decomposed_litter[:, 2] .> 0.0f0)

    invalid = init_soil(1, soilparams.soildepth, identity)
    invalid_crop = init_crop(1, identity)
    invalid_state = test_model_state(invalid_crop, invalid)
    invalid.thermal.temperature .= 10.0f0
    invalid.water.saturation_storage .= 100.0f0
    invalid.water.holding_capacity_storage .= 100.0f0
    invalid.water.relative_content .= 0.5f0
    invalid.carbon.fast .= -1.0f0
    invalid.carbon.slow .= -1.0f0
    invalid.nitrogen.fast .= -1.0f0
    invalid.nitrogen.slow .= -1.0f0

    soil_carbon!(invalid_state, invalid_state)
    soil_nitrogen!(invalid_state, invalid_state)

    # Match LPJmL's max(0, flux) guard: an invalid negative pool must not
    # create a reverse decomposition flux that increases respiration.
    @test all(iszero, invalid.carbon.decomposed_fast)
    @test all(iszero, invalid.carbon.decomposed_slow)
    @test all(iszero, invalid.nitrogen.decomposed_fast)
    @test all(iszero, invalid.nitrogen.decomposed_slow)
end
