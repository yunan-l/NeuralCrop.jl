using NeuralCrop
using Test

@testset "Crop nitrogen demand follows LPJmL ndemand_crop" begin
    crop = init_crop(1, identity)
    state = test_model_state(crop)
    crop.state.phenology.is_growing .= 1
    crop.state.carbon.leaf .= 2.0f0
    crop.state.carbon.root .= 3.0f0
    crop.state.carbon.pool .= 1.0f0
    crop.state.carbon.storage .= 4.0f0
    vcmax = Float32[10.0]
    temp = Float32[25.0]

    ndemand_crop!(state, cft1, vcmax, temp)

    expected_leaf = lpjmlparams.p * 1.0f-3 * vcmax[1] /
                    (86400.0f0 * 12.0f0 * 1.0f-6) +
                    cft1.ncleaf.low * crop.state.carbon.leaf[1]
    expected_nc = clamp(
        expected_leaf / crop.state.carbon.leaf[1],
        cft1.ncleaf.low,
        cft1.ncleaf.high,
    )
    expected_total = expected_leaf + expected_nc * (
        crop.state.carbon.root[1] / cft1.ratio.root +
        crop.state.carbon.pool[1] / cft1.ratio.pool +
        crop.state.carbon.storage[1] / cft1.ratio.sto
    )

    @test crop.auxiliary.stress.nitrogen_demand_leaf[1] ≈ expected_leaf atol = 1.0f-6
    @test crop.auxiliary.stress.nitrogen_demand_total[1] ≈ expected_total atol = 1.0f-6
    @test crop.auxiliary.stress.nitrogen_demand_total[1] >= crop.auxiliary.stress.nitrogen_demand_leaf[1]
end

@testset "Inactive crop has no nitrogen demand" begin
    crop = init_crop(1, identity)
    state = test_model_state(crop)
    crop.auxiliary.stress.nitrogen_demand_leaf .= 9.0f0
    crop.auxiliary.stress.nitrogen_demand_total .= 9.0f0

    ndemand_crop!(state, cft1, Float32[10.0], Float32[25.0])

    @test crop.auxiliary.stress.nitrogen_demand_leaf[1] == 0.0f0
    @test crop.auxiliary.stress.nitrogen_demand_total[1] == 0.0f0
end
