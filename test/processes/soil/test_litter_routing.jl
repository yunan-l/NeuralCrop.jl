using NeuralCrop
using Test

@testset "LPJmL litter spatial routing" begin
    @testset "sowing tillage conserves C and N" begin
        soil = init_soil(2, soilparams.soildepth, identity)
        crop = init_crop(2, identity)
        state = test_model_state(crop, soil)
        soil.management.tillage_fraction .= Float32[
            0.05 0 0
            0.95 1 0
            0 0 1
        ]
        soil.carbon.litter .= Float32[10 4; 2 3; 5 6]
        soil.nitrogen.litter .= Float32[1 0.4; 0.2 0.3; 0.5 0.6]
        crop.events.sowing .= Int32[1, 0]

        carbon_before = vec(sum(soil.carbon.litter, dims = 1))
        nitrogen_before = vec(sum(soil.nitrogen.litter, dims = 1))
        litter_tillage!(state, state)

        @test soil.carbon.litter[:, 1] ≈ Float32[0.5, 11.5, 5]
        @test soil.carbon.litter[:, 2] == Float32[4, 3, 6]
        @test soil.nitrogen.litter[:, 1] ≈ Float32[0.05, 1.15, 0.5]
        @test vec(sum(soil.carbon.litter, dims = 1)) ≈ carbon_before
        @test vec(sum(soil.nitrogen.litter, dims = 1)) ≈ nitrogen_before
        @test soil.management.tillage_carbon ≈ Float32[9.5, 0]
        @test soil.management.tillage_nitrogen ≈ Float32[0.95, 0]
    end

    @testset "daily bioturbation matches annual LPJmL fraction" begin
        soil = init_soil(1, soilparams.soildepth, identity)
        state = test_model_state(soil)
        soil.carbon.litter[:, 1] .= Float32[10, 2, 5]
        soil.nitrogen.litter[:, 1] .= Float32[1, 0.2, 0.5]
        carbon_before = sum(soil.carbon.litter)
        nitrogen_before = sum(soil.nitrogen.litter)

        expected_carbon_transfer = 10.0f0 * lpjmlparams.bioturbate
        expected_nitrogen_transfer = 1.0f0 * lpjmlparams.bioturbate
        litter_bioturbation!(state)

        @test soil.management.bioturbation_carbon[1] ≈ expected_carbon_transfer
        @test soil.management.bioturbation_nitrogen[1] ≈ expected_nitrogen_transfer
        @test sum(soil.carbon.litter) ≈ carbon_before atol = 2.0f-6
        @test sum(soil.nitrogen.litter) ≈ nitrogen_before atol = 2.0f-7
        @test soil.carbon.litter[3, 1] == 5.0f0
        @test (1.0f0 - lpjmlparams.bioturbate)^365 ≈ 0.5f0 atol = 2.0f-5
    end

    @testset "post-harvest setaside tills shoot but not root residues" begin
        soil = init_soil(1, soilparams.soildepth, identity)
        crop = init_crop(1, identity)
        output = init_output(1, identity)
        state = test_model_state(crop, soil; output)
        crop.state.carbon.leaf .= 2.0f0
        crop.state.carbon.pool .= 1.0f0
        crop.state.carbon.root .= 4.0f0
        crop.state.nitrogen.leaf .= 0.2f0
        crop.state.nitrogen.pool .= 0.1f0
        crop.state.nitrogen.root .= 0.4f0
        crop.state.phenology.harvesting_previous .= false
        crop.state.phenology.harvesting .= true

        harvest_crop!(state, state, output, Float32[0.5], 100)
        soil.management.tillage_fraction .= Float32[
            0.05 0 0
            0.95 1 0
            0 0 1
        ]
        NeuralCrop.route_harvest_carbon_input!(state, state)
        NeuralCrop.route_harvest_nitrogen_input!(state, state)

        @test soil.carbon.litter[:, 1] ≈ Float32[0.075, 1.425, 4]
        @test soil.nitrogen.litter[:, 1] ≈ Float32[0.0075, 0.1425, 0.4]
        @test soil.management.tillage_carbon[1] ≈ 1.425f0
        @test soil.management.tillage_nitrogen[1] ≈ 0.1425f0
    end

end
