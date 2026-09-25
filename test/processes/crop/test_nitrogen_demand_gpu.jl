using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA crop nitrogen demand" begin
    crop = init_crop(2, CuArray)
    state = test_model_state(crop)
    crop.state.phenology.is_growing .= 1
    crop.state.carbon.leaf .= 2.0f0
    crop.state.carbon.root .= 3.0f0
    crop.state.carbon.pool .= 1.0f0
    crop.state.carbon.storage .= 4.0f0

    ndemand_crop!(
        state,
        cft1,
        CuArray(Float32[10.0, 5.0]),
        CuArray(Float32[25.0, 15.0]),
    )

    leaf_demand = Array(crop.auxiliary.stress.nitrogen_demand_leaf)
    total_demand = Array(crop.auxiliary.stress.nitrogen_demand_total)
    @test all(isfinite, leaf_demand)
    @test all(total_demand .>= leaf_demand)
end
