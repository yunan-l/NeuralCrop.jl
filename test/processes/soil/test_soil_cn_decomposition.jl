using NeuralCrop
using Test

function initialize_soil_cn_case(litter_cn::T; mineral_n::T = T(2)) where {T <: AbstractFloat}
    soil = init_soil(T, 1, T.(soilparams.soildepth), identity)

    soil.carbon.litter .= reshape(T[30, 20, 10], 3, 1)
    soil.nitrogen.litter .= soil.carbon.litter ./ litter_cn
    soil.carbon.fast .= T(12)
    soil.carbon.slow .= T(24)
    soil.nitrogen.fast .= T(1.2)
    soil.nitrogen.slow .= T(1.2)
    soil.nitrogen.ammonium .= mineral_n / T(10)
    soil.nitrogen.nitrate .= mineral_n / T(10)

    soil.carbon.litter_response .= T(0.08)
    soil.nitrogen.litter_response .= T(0.08)
    soil.decomposition.shift_fast .= zero(T)
    soil.decomposition.shift_slow .= zero(T)
    soil.decomposition.shift_fast[1, 1] = one(T)
    soil.decomposition.shift_slow[1, 1] = one(T)

    soil.thermal.temperature .= T(10)
    soil.surface_litter.temperature .= T(10)
    soil.water.saturation_storage .= T(100)
    soil.water.holding_capacity_storage .= T(60)
    soil.water.wilting_storage .= T(10)
    soil.water.relative_content .= T(0.5)
    soil.properties.ph .= T(6.5)
    return soil
end

organic_nitrogen(soil) =
    sum(soil.nitrogen.litter) + sum(soil.nitrogen.fast) + sum(soil.nitrogen.slow)
mineral_nitrogen(soil) = sum(soil.nitrogen.ammonium) + sum(soil.nitrogen.nitrate)

@testset "Coupled LPJmL-style soil C-N decomposition" begin
    @testset "C and N use identical decay fractions" begin
        soil = initialize_soil_cn_case(10.0f0)
        state = test_model_state(soil)
        litter_carbon_before = copy(soil.carbon.litter)
        litter_nitrogen_before = copy(soil.nitrogen.litter)
        fast_carbon_before = copy(soil.carbon.fast)
        fast_nitrogen_before = copy(soil.nitrogen.fast)
        slow_carbon_before = copy(soil.carbon.slow)
        slow_nitrogen_before = copy(soil.nitrogen.slow)
        carbon_before = sum(soil.carbon.litter) + sum(soil.carbon.fast) + sum(soil.carbon.slow)
        nitrogen_before = organic_nitrogen(soil) + mineral_nitrogen(soil)

        # The coupled path uses one shared c_shift configuration for C and N.
        soil.nitrogen.litter_response .= 0.5f0
        @test !hasproperty(soil.carbon, :shift_fast)
        @test !hasproperty(soil.nitrogen, :shift_fast)

        soil_cn_decomposition!(state)

        @test soil.carbon.decomposed_litter ./ litter_carbon_before ≈
              soil.nitrogen.decomposed_litter ./ litter_nitrogen_before
        @test soil.carbon.decomposed_fast ./ fast_carbon_before ≈
              soil.nitrogen.decomposed_fast ./ fast_nitrogen_before
        @test soil.carbon.decomposed_slow ./ slow_carbon_before ≈
              soil.nitrogen.decomposed_slow ./ slow_nitrogen_before
        @test soil.nitrogen.litter_to_fast[1, 1] > 0.0f0
        @test soil.nitrogen.litter_to_fast[2, 1] == 0.0f0

        carbon_after = sum(soil.carbon.litter) + sum(soil.carbon.fast) + sum(soil.carbon.slow)
        nitrogen_after = organic_nitrogen(soil) + mineral_nitrogen(soil)
        @test carbon_after + sum(soil.carbon.heterotrophic_respiration) ≈
              carbon_before atol = 2.0f-5
        @test nitrogen_after + sum(soil.nitrogen.n2o_nitrification) ≈
              nitrogen_before atol = 2.0f-5
    end

    @testset "Litter C:N controls mineralization and immobilization" begin
        nitrogen_rich = initialize_soil_cn_case(5.0f0)
        soil_cn_decomposition!(test_model_state(nitrogen_rich))
        @test sum(nitrogen_rich.nitrogen.mineralization) > 0.0f0
        @test sum(nitrogen_rich.nitrogen.immobilization) == 0.0f0

        nitrogen_poor = initialize_soil_cn_case(80.0f0; mineral_n = 5.0f0)
        soil_cn_decomposition!(test_model_state(nitrogen_poor))
        @test sum(nitrogen_poor.nitrogen.immobilization) > 0.0f0

        supply_limited = initialize_soil_cn_case(80.0f0; mineral_n = 0.0f0)
        soil_cn_decomposition!(test_model_state(supply_limited))
        @test sum(supply_limited.nitrogen.immobilization) >= 0.0f0
        @test sum(supply_limited.nitrogen.immobilization) <
              sum(nitrogen_poor.nitrogen.immobilization)
        @test all(supply_limited.nitrogen.ammonium .>= 0.0f0)
        @test all(supply_limited.nitrogen.nitrate .>= 0.0f0)
    end

    @testset "Pre-crop and post-crop nitrogen stages are separated" begin
        soil = initialize_soil_cn_case(5.0f0; mineral_n = 0.0f0)
        state = test_model_state(soil)
        mineral_before = mineral_nitrogen(soil)
        soil_cn_decomposition!(state)
        mineral_after_decomposition = mineral_nitrogen(soil)

        @test mineral_after_decomposition > mineral_before
        @test all(iszero, soil.nitrogen.denitrification)
        @test all(iszero, soil.nitrogen.volatilization)

        # Represent plant uptake between LPJmL's pre- and post-crop stages.
        soil.nitrogen.ammonium .*= 0.5f0
        soil.nitrogen.nitrate .*= 0.5f0
        post_crop_nitrogen_losses!(
            state; air_temperature = Float32[20], wind_speed = Float32[2],
        )
        @test all(soil.nitrogen.denitrification .>= 0.0f0)
        @test all(soil.nitrogen.volatilization .>= 0.0f0)
        @test mineral_nitrogen(soil) <= mineral_after_decomposition
    end
end
