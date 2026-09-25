using NeuralCrop
using Test

@testset "LPJmL crop nitrogen limitation of Vcmax" begin
    crop = init_crop(2, identity)
    state = test_model_state(crop)
    crop.state.phenology.is_growing .= true
    crop.state.carbon.leaf .= 100.0f0
    crop.auxiliary.photosynthesis.potential_vcmax .= Float32[10, 10]
    crop.auxiliary.photosynthesis.vcmax .= crop.auxiliary.photosynthesis.potential_vcmax
    temperature = Float32[25, 25]

    # The demand equation and Vcmax-limit equation are analytical inverses
    # when optimal leaf N is available.
    ndemand_crop!(state, cft1, crop.auxiliary.photosynthesis.potential_vcmax, temperature)
    optimal_leaf_n = copy(crop.auxiliary.stress.nitrogen_demand_leaf)
    crop.auxiliary.stress.nitrogen_demand_leaf[2] = cft1.ncleaf.low * crop.state.carbon.leaf[2]
    limit_vcmax_by_nitrogen!(state, cft1, temperature)

    @test crop.auxiliary.photosynthesis.vcmax[1] ≈ 10.0f0 atol = 2.0f-5
    @test crop.auxiliary.photosynthesis.nitrogen_limitation[1] ≈ 1.0f0 atol = 2.0f-6
    @test 0.0f0 < crop.auxiliary.photosynthesis.vcmax[2] < crop.auxiliary.photosynthesis.potential_vcmax[2]
    @test crop.auxiliary.photosynthesis.nitrogen_limitation[2] < 1.0f-5
    @test crop.auxiliary.stress.nitrogen_demand_leaf[1] == optimal_leaf_n[1]
end
