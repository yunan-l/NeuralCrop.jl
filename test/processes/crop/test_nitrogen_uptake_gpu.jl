using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA separate NO3/NH4 uptake conserves nitrogen" begin
    crop = init_crop(1, CuArray)
    soil = init_soil(1, soilparams.soildepth, CuArray)
    state = test_model_state(crop, soil)

    crop.state.phenology.is_growing .= 1
    crop.state.nitrogen.total .= 0.1f0
    crop.state.carbon.leaf .= 20.0f0
    crop.state.carbon.root .= 100.0f0
    crop.auxiliary.stress.nitrogen_demand_leaf .= 0.4f0
    crop.auxiliary.stress.nitrogen_demand_total .= 1.0f0
    crop.auxiliary.root.distribution .= CuArray(Float32[0.5, 0.3, 0.2, 0.0, 0.0])
    soil.water.relative_content .= 1.0f0
    soil.water.saturation_fraction .= 0.4f0
    soil.thermal.temperature .= 15.0f0
    soil.nitrogen.nitrate .= CuArray(reshape(Float32[1.0, 0.0, 0.5, 0.0, 0.0], 5, 1))
    soil.nitrogen.ammonium .= CuArray(reshape(Float32[0.0, 0.7, 0.0, 0.0, 0.0], 5, 1))

    plant_before = Array(crop.state.nitrogen.total)[1]
    soil_before = sum(Array(soil.nitrogen.nitrate)) + sum(Array(soil.nitrogen.ammonium))
    nuptake_crop!(state, cft1, state)

    uptake = Array(crop.fluxes.nitrogen.uptake)[1]
    plant_gain = Array(crop.state.nitrogen.total)[1] - plant_before
    soil_loss = soil_before - sum(Array(soil.nitrogen.nitrate)) - sum(Array(soil.nitrogen.ammonium))

    @test uptake > 0.0f0
    @test plant_gain ≈ uptake atol = 1.0f-6
    @test soil_loss ≈ uptake atol = 1.0f-6
    @test all(Array(soil.nitrogen.nitrate) .>= 0.0f0)
    @test all(Array(soil.nitrogen.ammonium) .>= 0.0f0)
end

@testset "CUDA automatic fertilizer supplies the remaining demand" begin
    crop = init_crop(1, CuArray)
    soil = init_soil(1, soilparams.soildepth, CuArray)
    state = test_model_state(crop, soil)
    crop.state.phenology.is_growing .= 1
    crop.state.nitrogen.total .= 0.1f0
    crop.state.carbon.leaf .= 20.0f0
    crop.state.carbon.root .= 100.0f0
    crop.auxiliary.stress.nitrogen_demand_leaf .= 0.4f0
    crop.auxiliary.stress.nitrogen_demand_total .= 1.0f0

    nuptake_crop!(state, cft1, state; auto_fertilizer = true)

    @test Array(crop.state.nitrogen.total)[1] ≈ 1.0f0 atol = 1.0f-6
    @test Array(crop.fluxes.nitrogen.uptake)[1] ≈ 0.9f0 atol = 1.0f-6
    @test Array(crop.fluxes.nitrogen.auto_fertilizer)[1] ≈ 0.9f0 atol = 1.0f-6
    @test Array(crop.state.nitrogen.sufficiency)[1] == 1.0f0
end
