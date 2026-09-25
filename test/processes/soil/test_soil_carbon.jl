using NeuralCrop
using Test

@testset "Conservative soil-carbon decomposition" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)
    state = test_model_state(crop, soil)

    soil.carbon.litter .= 2.0f0
    soil.carbon.fast .= 3.0f0
    soil.carbon.slow .= 4.0f0
    soil.carbon.litter_response .= Float32[0.97, 0.97, 0.30]
    soil.decomposition.shift_fast[1, 1] = 1.0f0
    soil.decomposition.shift_slow[1, 1] = 1.0f0
    soil.thermal.temperature .= 10.0f0
    soil.water.relative_content .= 0.5f0
    soil.water.holding_capacity_storage .= 80.0f0
    soil.water.wilting_storage .= 10.0f0
    soil.water.saturation_storage .= 100.0f0

    carbon_before = sum(soil.carbon.litter) +
                    sum(soil.carbon.fast) + sum(soil.carbon.slow)
    soil_carbon!(state, state)
    carbon_after = sum(soil.carbon.litter) +
                   sum(soil.carbon.fast) + sum(soil.carbon.slow)

    @test carbon_after < carbon_before
    @test carbon_after + sum(soil.carbon.heterotrophic_respiration) ≈
          carbon_before atol = 1.0f-5
    @test soil.carbon.heterotrophic_respiration[1] ≈
          sum(soil.carbon.decomposed_litter) * lpjmlparams.atmfrac +
          sum(soil.carbon.decomposed_fast) + sum(soil.carbon.decomposed_slow)
    @test all(soil.carbon.litter .>= 0.0f0)
    @test all(soil.carbon.fast .>= 0.0f0)
    @test all(soil.carbon.slow .>= 0.0f0)
end

@testset "Daily carbon-balance diagnostics" begin
    crop = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    state = test_model_state(crop, soil)
    balance = @inferred init_carbon_balance(4, 1, identity)

    crop.state.carbon.leaf .= 1.0f0
    soil.carbon.litter .= 2.0f0
    soil.carbon.fast .= 3.0f0
    soil.carbon.slow .= 4.0f0
    NeuralCrop.record_carbon_balance_start!(balance, 1, state, state)
    NeuralCrop.record_carbon_balance_end!(balance, 1, state, state)

    @test balance.plant_before[1, 1] == 1.0f0
    @test balance.soil_before[1, 1] == 41.0f0
    @test balance.residual[1, 1] ≈ 0.0f0 atol = 1.0f-6

    # NPP is the net atmospheric input to the tracked plant-plus-soil system.
    NeuralCrop.record_carbon_balance_start!(balance, 2, state, state)
    crop.state.carbon.leaf .+= 2.0f0
    crop.fluxes.carbon.npp .= 2.0f0
    NeuralCrop.record_carbon_balance_end!(balance, 2, state, state)
    @test balance.net_primary_production[2, 1] == 2.0f0
    @test balance.residual[2, 1] ≈ 0.0f0 atol = 1.0f-6

    # Harvest residue transfer is internal; only harvested material leaves.
    crop.state.carbon.leaf .= 0.2f0
    crop.state.carbon.root .= 0.3f0
    crop.state.carbon.storage .= 0.3f0
    crop.state.carbon.pool .= 0.2f0
    crop.fluxes.carbon.npp .= 0.0f0
    crop.state.phenology.harvesting_previous .= false
    crop.state.phenology.harvesting .= true
    soil.carbon.litter .= 0.0f0
    soil.carbon.fast .= 0.0f0
    soil.carbon.slow .= 0.0f0
    soil.management.tillage_fraction .= Float32[1 0 0; 0 1 0; 0 0 1]
    output = init_output(1, identity)

    NeuralCrop.record_carbon_balance_start!(balance, 3, state, state)
    harvest_crop!(state, state, output, Float32[0.5], 100)
    NeuralCrop.record_carbon_balance_after_harvest!(
        balance, 3, state, state, Float32[0.5],
    )
    soil_carbon!(state, state)
    crop.state.carbon.leaf .= 0.0f0
    crop.state.carbon.root .= 0.0f0
    crop.state.carbon.storage .= 0.0f0
    crop.state.carbon.pool .= 0.0f0
    NeuralCrop.record_carbon_balance_end!(balance, 3, state, state)

    @test balance.residue_transfer[3, 1] ≈ 0.5f0 atol = 1.0f-6
    @test balance.harvest_export[3, 1] ≈ 0.5f0 atol = 1.0f-6
    @test balance.residual[3, 1] ≈ 0.0f0 atol = 1.0f-6

    # Manure carbon is an external input, not unexplained soil-C creation.
    NeuralCrop.record_carbon_balance_start!(balance, 4, state, state)
    crop.fluxes.nitrogen.prescribed_manure_input .= 1.0f0
    soil.carbon.litter[2, 1] += lpjmlparams.manure_cn
    NeuralCrop.record_carbon_balance_after_cultivate!(balance, 4, state)
    NeuralCrop.record_carbon_balance_end!(balance, 4, state, state)
    @test balance.manure_input[4, 1] == lpjmlparams.manure_cn
    @test balance.residual[4, 1] ≈ 0.0f0 atol = 1.0f-6
end
