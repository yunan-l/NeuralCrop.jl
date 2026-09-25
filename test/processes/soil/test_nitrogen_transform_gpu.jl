using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA soil nitrogen transformations" begin
    cell_size = 2
    soil = init_soil(cell_size, soilparams.soildepth, CuArray)
    state = test_model_state(soil)
    soil.properties.ph .= 6.5f0
    soil.nitrogen.ammonium .= 0.2f0
    soil.nitrogen.nitrate .= 0.1f0
    soil.nitrogen.fast .= 2.0f0
    soil.nitrogen.slow .= 4.0f0
    soil.decomposition.shift_fast[1, :] .= 1.0f0
    soil.decomposition.shift_slow[1, :] .= 1.0f0
    soil.water.relative_content .= 0.7f0
    soil.water.holding_capacity_storage .= 100.0f0
    soil.water.saturation_storage .= 120.0f0
    soil.thermal.temperature .= 20.0f0

    stock_before = sum(Array(soil.nitrogen.ammonium .+ soil.nitrogen.nitrate .+
                             soil.nitrogen.fast .+ soil.nitrogen.slow))
    # `nitrogen_transform!` receives SOM decomposition fluxes after the source
    # pools have already been reduced by `soil_nitrogen!`.
    soil.nitrogen.decomposed_fast .= 0.001f0
    soil.nitrogen.decomposed_slow .= 0.001f0
    soil.nitrogen.fast .-= soil.nitrogen.decomposed_fast
    soil.nitrogen.slow .-= soil.nitrogen.decomposed_slow
    nitrogen_transform!(
        state;
        air_temperature = CuArray(Float32[20, 21]),
        wind_speed = CuArray(Float32[1.5, 2.0]),
    )
    stock_after = sum(Array(soil.nitrogen.ammonium .+ soil.nitrogen.nitrate .+
                            soil.nitrogen.fast .+ soil.nitrogen.slow))
    gases = sum(Array(soil.nitrogen.n2o_nitrification .+
                       soil.nitrogen.n2o_denitrification .+
                       soil.nitrogen.n2_denitrification)) +
            sum(Array(soil.nitrogen.volatilization))

    @test isapprox(stock_before, stock_after + gases; atol = 2.0f-5)
    @test all(isfinite, Array(soil.nitrogen.ammonium))
    @test all(Array(soil.nitrogen.ammonium) .>= 0.0f0)
    @test all(Array(soil.nitrogen.nitrate) .>= 0.0f0)
end
