using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("CUDA is required for this independent GPU test")

@testset "GPU litter spatial routing" begin
    soil = init_soil(2, CuArray(soilparams.soildepth), CuArray)
    crop = init_crop(2, CuArray)
    state = test_model_state(crop, soil)
    soil.management.tillage_fraction .= CuArray(Float32[
        0.05 0 0
        0.95 1 0
        0 0 1
    ])
    soil.carbon.litter .= CuArray(Float32[10 4; 2 3; 5 6])
    soil.nitrogen.litter .= CuArray(Float32[1 0.4; 0.2 0.3; 0.5 0.6])
    crop.events.sowing .= CuArray(Int32[1, 0])

    carbon_before = Array(sum(soil.carbon.litter, dims = 1))
    nitrogen_before = Array(sum(soil.nitrogen.litter, dims = 1))
    litter_tillage!(state, state)
    litter_bioturbation!(state)

    @test Array(sum(soil.carbon.litter, dims = 1)) ≈ carbon_before atol = 2.0f-6
    @test Array(sum(soil.nitrogen.litter, dims = 1)) ≈ nitrogen_before atol = 2.0f-7
    @test all(Array(soil.management.tillage_carbon) .>= 0.0f0)
    @test all(Array(soil.management.bioturbation_carbon) .>= 0.0f0)
    @test Array(soil.carbon.litter[3:3, :]) == Float32[5 6]
end
