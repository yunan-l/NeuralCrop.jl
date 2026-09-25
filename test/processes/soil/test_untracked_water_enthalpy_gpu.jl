using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA LPJmL untracked water-mass enthalpy" begin
    cells = 32
    soil = init_soil(cells, soilparams.soildepth, CuArray)
    crop = init_crop(cells, CuArray)
    state = test_model_state(crop, soil)
    soil.properties.sand_fraction .= 0.4f0
    soil.properties.clay_fraction .= 0.2f0
    soil.water.storage .= CuArray(repeat(Float32[40, 60, 100, 200, 200], 1, cells))
    pedotransfer!(state)
    soil_temperature!(state, CUDA.fill(5.0f0, cells), CUDA.fill(5.0f0, cells))

    crop.fluxes.water.transpiration_layer .= 0.5f0
    soil.water.evaporation .= 0.25f0
    soil_evapotranspiration!(state, state)
    soil_temperature!(state, CUDA.fill(5.0f0, cells))
    synchronize()

    expected = -3.75f0 * 0.001f0 *
        (soil_thermal_params.volumetric_fusion_heat +
         soil_thermal_params.water_heat_capacity * 5.0f0)
    flux = Array(soil.thermal.untracked_water_energy_flux)
    temperature = Array(soil.thermal.temperature)
    @test all(isapprox.(flux, expected; rtol = 2.0f-5))
    @test maximum(abs.(temperature .- 5.0f0)) < 2.0f-5
    @test all(isfinite, Array(soil.thermal.enthalpy))
end
