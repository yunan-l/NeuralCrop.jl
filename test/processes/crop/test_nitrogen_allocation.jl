using NeuralCrop
using Test

@testset "Crop organ nitrogen is a redistributed stock" begin
    crop = init_crop(1, identity)
    state = test_model_state(crop)
    crop.state.phenology.is_growing .= 1
    crop.state.nitrogen.total .= 0.7f0
    crop.state.carbon.leaf .= 2.0f0
    crop.state.carbon.root .= 3.0f0
    crop.state.carbon.storage .= 4.0f0
    crop.state.carbon.pool .= 1.0f0

    # Deliberately inconsistent old organ values must not be accumulated.
    crop.state.nitrogen.leaf .= 10.0f0
    crop.state.nitrogen.root .= 20.0f0
    crop.state.nitrogen.storage .= 30.0f0
    crop.state.nitrogen.pool .= 40.0f0

    NeuralCrop.allocate_crop_nitrogen!(state, cft1)
    first_allocation = Float32[
        crop.state.nitrogen.leaf[1], crop.state.nitrogen.root[1], crop.state.nitrogen.storage[1], crop.state.nitrogen.pool[1],
    ]

    @test sum(first_allocation) ≈ crop.state.nitrogen.total[1] atol = 1.0f-6
    @test all(first_allocation .>= 0.0f0)

    NeuralCrop.allocate_crop_nitrogen!(state, cft1)
    second_allocation = Float32[
        crop.state.nitrogen.leaf[1], crop.state.nitrogen.root[1], crop.state.nitrogen.storage[1], crop.state.nitrogen.pool[1],
    ]

    @test second_allocation ≈ first_allocation atol = 1.0f-7
    @test sum(second_allocation) ≈ crop.state.nitrogen.total[1] atol = 1.0f-6
end

@testset "Zero leaf carbon clears organ nitrogen safely" begin
    crop = init_crop(1, identity)
    state = test_model_state(crop)
    crop.state.phenology.is_growing .= 1
    crop.state.nitrogen.total .= 0.7f0
    crop.state.carbon.leaf .= 0.0f0
    crop.state.carbon.root .= 3.0f0
    crop.state.nitrogen.leaf .= 1.0f0
    crop.state.nitrogen.root .= 1.0f0

    NeuralCrop.allocate_crop_nitrogen!(state, cft1)

    @test crop.state.nitrogen.leaf[1] == 0.0f0
    @test crop.state.nitrogen.root[1] == 0.0f0
    @test crop.state.nitrogen.storage[1] == 0.0f0
    @test crop.state.nitrogen.pool[1] == 0.0f0
end

@testset "Integrated automatic-fertilizer nitrogen cycle" begin
    crop = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    state = test_model_state(crop, soil)
    crop.state.phenology.is_growing .= 1
    crop.state.nitrogen.total .= 0.1f0
    crop.state.carbon.leaf .= 2.0f0
    crop.state.carbon.root .= 3.0f0
    crop.state.carbon.pool .= 1.0f0
    crop.state.carbon.storage .= 4.0f0

    crop_nitrogen!(
        state,
        cft1,
        state,
        Float32[10.0],
        Float32[25.0];
        auto_fertilizer = true,
    )

    organ_n = crop.state.nitrogen.leaf[1] + crop.state.nitrogen.root[1] + crop.state.nitrogen.pool[1] + crop.state.nitrogen.storage[1]
    @test crop.state.nitrogen.total[1] ≈ crop.auxiliary.stress.nitrogen_demand_total[1] atol = 1.0f-6
    @test organ_n ≈ crop.state.nitrogen.total[1] atol = 1.0f-6
    @test crop.fluxes.nitrogen.auto_fertilizer[1] > 0.0f0
    @test crop.state.nitrogen.sufficiency[1] == 1.0f0
end
