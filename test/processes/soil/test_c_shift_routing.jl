using NeuralCrop
using Test

@testset "Fixed post-spin-up c_shift routing" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)
    state = test_model_state(crop, soil)
    fast_shift = Float32[0.40, 0.25, 0.15, 0.10, 0.10]
    slow_shift = Float32[0.55, 0.20, 0.10, 0.10, 0.05]

    soil.decomposition.shift_fast[:, 1] .= fast_shift
    soil.decomposition.shift_slow[:, 1] .= slow_shift
    soil.carbon.litter[NeuralCrop.ROOT_LITTER, 1] = 10.0f0
    soil.nitrogen.litter[NeuralCrop.ROOT_LITTER, 1] = 0.4f0
    soil.carbon.litter_response .= 0.0f0
    soil.nitrogen.litter_response .= 0.0f0
    soil.carbon.litter_response[NeuralCrop.ROOT_LITTER] = 1.0f0
    soil.nitrogen.litter_response[NeuralCrop.ROOT_LITTER] = 1.0f0
    soil.thermal.temperature .= 10.0f0
    soil.water.relative_content .= 0.5f0
    soil.water.holding_capacity_storage .= 80.0f0
    soil.water.wilting_storage .= 10.0f0
    soil.water.saturation_storage .= 100.0f0
    soil.nitrogen.ammonium .= 0.2f0
    soil.nitrogen.nitrate .= 0.1f0

    @test sum(soil.decomposition.shift_fast; dims = 1)[1] ≈ 1.0f0
    @test sum(soil.decomposition.shift_slow; dims = 1)[1] ≈ 1.0f0

    soil_carbon!(state, state)
    soil_nitrogen!(state, state; air_temperature = Float32[10], wind_speed = Float32[1.5])

    decomposed_carbon = sum(soil.carbon.decomposed_litter)
    decomposed_nitrogen = sum(soil.nitrogen.decomposed_litter)
    retained_carbon = decomposed_carbon * (1.0f0 - lpjmlparams.atmfrac)
    retained_nitrogen = decomposed_nitrogen * (1.0f0 - lpjmlparams.atmfrac)

    @test soil.carbon.litter_to_fast[:, 1] ≈
          fast_shift .* retained_carbon .* lpjmlparams.fastfrac
    @test soil.carbon.litter_to_slow[:, 1] ≈
          slow_shift .* retained_carbon .* (1.0f0 - lpjmlparams.fastfrac)
    @test soil.nitrogen.litter_to_fast[:, 1] ≈
          fast_shift .* retained_nitrogen .* lpjmlparams.fastfrac
    @test soil.nitrogen.litter_to_slow[:, 1] ≈
          slow_shift .* retained_nitrogen .* (1.0f0 - lpjmlparams.fastfrac)
    @test sum(soil.carbon.litter_to_fast .+ soil.carbon.litter_to_slow) ≈
          retained_carbon atol = 1.0f-6
    @test sum(soil.nitrogen.litter_to_fast .+ soil.nitrogen.litter_to_slow) ≈
          retained_nitrogen atol = 1.0f-7
    @test sum(root_distribution(0.96f0)) ≈ 1.0 atol = 1.0e-6

    response_sum = reshape(Float32[1, 2, 3, 4, 5], 5, 1)
    calibrated_fast = similar(response_sum)
    calibrated_slow = similar(response_sum)
    equilibrated_c_shift!(
        calibrated_fast, calibrated_slow, response_sum,
        reshape(fast_shift, :, 1), reshape(slow_shift, :, 1),
    )
    @test sum(calibrated_fast; dims = 1)[1] ≈ 1.0f0
    @test sum(calibrated_slow; dims = 1)[1] ≈ 1.0f0
    @test calibrated_fast[:, 1] ≈ fast_shift .* response_sum[:, 1] ./
        sum(fast_shift .* response_sum[:, 1])
    @test calibrated_slow[:, 1] ≈ slow_shift .* response_sum[:, 1] ./
        sum(slow_shift .* response_sum[:, 1])
end
