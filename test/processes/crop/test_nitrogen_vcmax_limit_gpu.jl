using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA crop nitrogen limitation of Vcmax" begin
    crop = init_crop(2, CuArray)
    state = test_model_state(crop)
    crop.state.phenology.is_growing .= true
    crop.state.carbon.leaf .= 100.0f0
    crop.auxiliary.photosynthesis.potential_vcmax .= CuArray(Float32[10, 10])
    crop.auxiliary.photosynthesis.vcmax .= crop.auxiliary.photosynthesis.potential_vcmax
    temperature = CuArray(Float32[25, 25])

    ndemand_crop!(state, cft1, crop.auxiliary.photosynthesis.potential_vcmax, temperature)
    @views crop.auxiliary.stress.nitrogen_demand_leaf[2:2] .= cft1.ncleaf.low .* crop.state.carbon.leaf[2:2]
    limit_vcmax_by_nitrogen!(state, cft1, temperature)

    vcmax = Array(crop.auxiliary.photosynthesis.vcmax)
    limitation = Array(crop.auxiliary.photosynthesis.nitrogen_limitation)
    @test vcmax[1] ≈ 10.0f0 atol = 2.0f-5
    @test 0.0f0 < vcmax[2] < 10.0f0
    @test limitation[1] ≈ 1.0f0 atol = 2.0f-6
    @test limitation[2] < 1.0f-5
end
