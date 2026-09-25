using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA fixed post-spin-up c_shift routing" begin
    cells = 32
    soil = init_soil(cells, soilparams.soildepth, CuArray)
    crop = init_crop(cells, CuArray)
    state = test_model_state(crop, soil)
    fast_shift = Float32[0.40, 0.25, 0.15, 0.10, 0.10]
    slow_shift = Float32[0.55, 0.20, 0.10, 0.10, 0.05]

    soil.decomposition.shift_fast .= CuArray(repeat(fast_shift, 1, cells))
    soil.decomposition.shift_slow .= CuArray(repeat(slow_shift, 1, cells))
    root_litter = NeuralCrop.ROOT_LITTER
    soil.carbon.litter[root_litter:root_litter, :] .= 10.0f0
    soil.nitrogen.litter[root_litter:root_litter, :] .= 0.4f0
    soil.carbon.litter_response .= CuArray(Float32[0, 0, 1])
    soil.nitrogen.litter_response .= CuArray(Float32[0, 0, 1])
    soil.thermal.temperature .= 10.0f0
    soil.water.relative_content .= 0.5f0
    soil.water.holding_capacity_storage .= 80.0f0
    soil.water.wilting_storage .= 10.0f0
    soil.water.saturation_storage .= 100.0f0
    soil.nitrogen.ammonium .= 0.2f0
    soil.nitrogen.nitrate .= 0.1f0

    soil_carbon!(state, state)
    soil_nitrogen!(
        state, state;
        air_temperature = CUDA.fill(10.0f0, cells),
        wind_speed = CUDA.fill(1.5f0, cells),
    )
    synchronize()

    carbon_decomposed = Array(sum(soil.carbon.decomposed_litter; dims = 1))
    carbon_routed = Array(sum(soil.carbon.litter_to_fast .+
                              soil.carbon.litter_to_slow; dims = 1))
    nitrogen_decomposed = Array(sum(soil.nitrogen.decomposed_litter; dims = 1))
    nitrogen_routed = Array(sum(soil.nitrogen.litter_to_fast .+
                                soil.nitrogen.litter_to_slow; dims = 1))

    @test carbon_routed ≈ carbon_decomposed .* (1.0f0 - lpjmlparams.atmfrac) atol = 2.0f-6
    @test nitrogen_routed ≈ nitrogen_decomposed .* (1.0f0 - lpjmlparams.atmfrac) atol = 2.0f-7
    @test all(isfinite, Array(soil.carbon.fast))
    @test all(isfinite, Array(soil.nitrogen.fast))
end
