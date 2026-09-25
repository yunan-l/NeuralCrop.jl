using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA five-layer implicit soil heat conduction" begin
    cell_size = 32
    soil = init_soil(cell_size, soilparams.soildepth, CuArray)
    state = test_model_state(soil)
    soil.thermal.diffusivity_0 .= 0.6f0
    soil.thermal.diffusivity_15 .= 0.7f0
    soil.water.relative_content .= 0.1f0
    litter_carbon_2cm = 20.0f0 * soil_thermal_params.litter_carbon_fraction *
                       soil_thermal_params.litter_bulk_density
    soil.carbon.litter[1, 9:16] .= litter_carbon_2cm
    soil.snow.height[17:24] .= 0.67f0
    soil.carbon.litter[1, 25:32] .= litter_carbon_2cm
    soil.snow.height[25:32] .= 0.67f0
    update_surface_litter_properties!(state)

    air_temperature = CuArray(fill(30.0f0, cell_size))
    initial_temperature = CuArray(fill(10.0f0, cell_size))
    soil_temperature!(state, air_temperature, initial_temperature)
    synchronize()

    profile = Array(soil.thermal.temperature)
    @test all(Array(soil.thermal.initialized))
    @test all(isfinite, profile)
    @test all(profile[1, :] .> profile[2, :])
    @test all(profile[2, :] .> profile[3, :])
    @test all(profile[3, :] .> profile[4, :])
    @test all(profile[4, :] .> profile[5, :])
    @test all((profile .>= 10.0f0) .& (profile .<= 30.0f0))
    @test all(profile[1, 9:16] .< profile[1, 1:8])
    @test all(profile[1, 17:24] .< profile[1, 1:8])
    @test all(profile[1, 25:32] .< profile[1, 9:16])
    @test all(profile[1, 25:32] .< profile[1, 17:24])
end
