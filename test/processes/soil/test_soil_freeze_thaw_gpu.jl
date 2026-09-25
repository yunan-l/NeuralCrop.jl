using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA five-layer freeze-thaw" begin
    cells = 32
    soil = init_soil(cells, soilparams.soildepth, CuArray)
    state = test_model_state(soil)
    soil.water.storage .= 50.0f0
    pedotransfer!(state)
    total_before = Array(vec(sum(
        soil.water.storage + soil.water.ice_storage; dims = 1,
    )))

    soil_temperature!(state, CUDA.fill(-20.0f0, cells), CUDA.fill(2.0f0, cells))
    total_after = Array(vec(sum(
        soil.water.storage + soil.water.ice_storage; dims = 1,
    )))
    frozen_fraction = Array(soil.thermal.frozen_fraction)
    @test total_after ≈ total_before atol = 2.0f-5
    @test all(isfinite, Array(soil.thermal.enthalpy))
    @test all((0.0f0 .<= frozen_fraction) .& (frozen_fraction .<= 1.0f0))
    @test all(Array(vec(sum(soil.water.ice_storage; dims = 1))) .> 0.0f0)
    relative_residual = abs.(Array(soil.thermal.energy_residual)) ./
        max.(abs.(Array(soil.thermal.surface_energy_flux)), 1.0f0)
    @test maximum(relative_residual) < 2.0f-5
end
