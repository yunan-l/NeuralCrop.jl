using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA crop organ nitrogen redistribution" begin
    crop = init_crop(2, CuArray)
    state = test_model_state(crop)
    crop.state.phenology.is_growing .= 1
    crop.state.nitrogen.total .= 0.7f0
    crop.state.carbon.leaf .= 2.0f0
    crop.state.carbon.root .= 3.0f0
    crop.state.carbon.storage .= 4.0f0
    crop.state.carbon.pool .= 1.0f0
    crop.state.nitrogen.leaf .= 10.0f0
    crop.state.nitrogen.root .= 20.0f0
    crop.state.nitrogen.storage .= 30.0f0
    crop.state.nitrogen.pool .= 40.0f0

    NeuralCrop.allocate_crop_nitrogen!(state, cft1)
    first_sum = Array(crop.state.nitrogen.leaf .+ crop.state.nitrogen.root .+ crop.state.nitrogen.storage .+ crop.state.nitrogen.pool)
    first_leafn = Array(crop.state.nitrogen.leaf)

    NeuralCrop.allocate_crop_nitrogen!(state, cft1)
    second_sum = Array(crop.state.nitrogen.leaf .+ crop.state.nitrogen.root .+ crop.state.nitrogen.storage .+ crop.state.nitrogen.pool)
    second_leafn = Array(crop.state.nitrogen.leaf)

    @test first_sum ≈ fill(0.7f0, 2) atol = 1.0f-6
    @test second_sum ≈ first_sum atol = 1.0f-6
    @test second_leafn ≈ first_leafn atol = 1.0f-7
end
