using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA LPJmL daily soil decomposition flux guards" begin
    cells = 32
    crop = init_crop(cells, CuArray)
    soil = init_soil(cells, soilparams.soildepth, CuArray)
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
    soil.thermal.temperature .= -20.0f0

    soil_carbon!(state, state)
    soil_nitrogen!(state, state)
    synchronize()

    @test all(iszero, Array(soil.carbon.decomposed_litter))
    @test all(iszero, Array(soil.nitrogen.decomposed_litter))
    @test all(iszero, Array(soil.carbon.decomposed_fast))
    @test all(iszero, Array(soil.carbon.decomposed_slow))
    @test all(iszero, Array(soil.nitrogen.decomposed_fast))
    @test all(iszero, Array(soil.nitrogen.decomposed_slow))

    soil.thermal.temperature .= 10.0f0
    soil.carbon.litter .= 0.0f0
    soil.nitrogen.litter .= 0.0f0
    soil.carbon.fast .= -1.0f0
    soil.carbon.slow .= -1.0f0
    soil.nitrogen.fast .= -1.0f0
    soil.nitrogen.slow .= -1.0f0
    soil_carbon!(state, state)
    soil_nitrogen!(state, state)
    synchronize()

    @test all(iszero, Array(soil.carbon.decomposed_fast))
    @test all(iszero, Array(soil.carbon.decomposed_slow))
    @test all(iszero, Array(soil.nitrogen.decomposed_fast))
    @test all(iszero, Array(soil.nitrogen.decomposed_slow))
end
