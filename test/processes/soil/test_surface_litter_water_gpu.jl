using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA surface-litter hydrology" begin
    cells = 32
    soil = init_soil(cells, soilparams.soildepth, CuArray)
    crop = init_crop(cells, CuArray)
    state = test_model_state(crop, soil)
    soil.carbon.litter[1, :] .= 20.0f0 * 0.42f0 * 71.1f0
    update_surface_litter_properties!(state)

    soil.water.infiltration .= 10.0f0
    surface_litter_interception!(state)
    @test all(Array(soil.surface_litter.water_storage) .≈
              Array(soil.surface_litter.water_capacity))
    @test all(Array(soil.water.infiltration) .>= 0.0f0)

    soil.water.storage .= 50.0f0
    soil.water.wilting_storage .= 10.0f0
    soil.water.holding_capacity_storage .= 100.0f0
    evaporation!(CUDA.fill(2.0f0, cells), state, state)
    @test all(isfinite, Array(soil.surface_litter.evaporation))
    @test all(Array(soil.surface_litter.evaporation) .> 0.0f0)
    @test all(Array(soil.surface_litter.water_storage) .>= 0.0f0)
end
