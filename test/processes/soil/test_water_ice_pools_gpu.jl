using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA LPJmL three-reservoir water-ice partition" begin
    cells = 32
    soil = init_soil(cells, soilparams.soildepth, CuArray)
    state = test_model_state(soil)
    soil.water.storage .= 50.0f0
    pedotransfer!(state)
    soil_temperature!(
        state,
        CUDA.fill(-20.0f0, cells),
        CUDA.fill(2.0f0, cells),
    )
    diagnostics = init_thermal_balance(1, cells, CuArray)
    NeuralCrop.record_thermal_balance!(diagnostics, 1, state)

    component_ice = Array(
        soil.water.wilting_ice_fraction .* soil.water.wilting_storage .+
        soil.water.available_ice_storage .+
        soil.water.free_ice_storage
    )
    total_ice = Array(soil.water.ice_storage)
    @test component_ice ≈ total_ice atol = 2.0f-5
    @test all(isfinite, component_ice)
    @test all(Array(soil.water.free_ice_storage)[
        Array(soil.water.wilting_ice_fraction) .< 1.0f0
    ] .== 0.0f0)
    @test maximum(Array(diagnostics.ice_pool_residual)) <= 2.0f-5
end
