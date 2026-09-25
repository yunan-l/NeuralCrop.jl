using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA LPJmL percolation enthalpy" begin
    cells = 32
    soil = init_soil(cells, soilparams.soildepth, CuArray)
    crop = init_crop(cells, CuArray)
    state = test_model_state(crop, soil)
    soil.properties.sand_fraction .= 0.4f0
    soil.properties.clay_fraction .= 0.2f0
    soil.water.storage .= CuArray(repeat(Float32[40, 60, 100, 200, 200], 1, cells))
    crop.fluxes.water.interception .= 0.0f0
    pedotransfer!(state)
    soil_temperature!(state, CUDA.fill(10.0f0, cells), CUDA.fill(10.0f0, cells))
    soil_infiltration!(
        state, state, CUDA.fill(4.0f0, cells);
        snowmelt = CUDA.fill(2.0f0, cells),
        air_temperature = CUDA.fill(10.0f0, cells),
    )
    synchronize()

    rain_energy = Array(soil.thermal.rain_energy_input)
    melt_energy = Array(soil.thermal.snowmelt_energy_input)
    residual = Array(soil.thermal.percolation_energy_residual)
    temperatures = Array(soil.thermal.temperature)
    @test all(rain_energy .> 0.0f0)
    @test all(melt_energy .> 0.0f0)
    @test all(temperatures[1, :] .< 10.0f0)
    @test all(isfinite, temperatures)
    @test maximum(abs, residual) < 2.0f0
    @test all(iszero, Array(soil.thermal.percolation_energy))
end
