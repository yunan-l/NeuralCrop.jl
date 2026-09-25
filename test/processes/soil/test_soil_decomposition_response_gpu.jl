using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA LPJmL soil decomposition response" begin
    cells = 32
    soil = init_soil(cells, soilparams.soildepth, CuArray)
    state = test_model_state(soil)
    soil.thermal.temperature .= 10.0f0
    soil.water.saturation_storage .= 100.0f0
    soil.water.holding_capacity_storage .= 100.0f0
    soil.water.relative_content .= 0.5f0
    soil.surface_litter.temperature .= -20.0f0
    soil.surface_litter.water_capacity .= 1.0f0
    soil.surface_litter.water_storage .= 0.5f0
    soil_decomp_response!(state)
    synchronize()

    response = Array(soil.decomposition.response)
    litter_response = Array(soil.decomposition.litter_response)
    @test all(isfinite, response)
    @test all((response .>= 0.0f0) .& (response .<= 1.0f0))
    @test all(iszero, litter_response[1, :])
    @test litter_response[2, :] ≈ response[1, :]
    @test litter_response[3, :] ≈ response[1, :]
end
