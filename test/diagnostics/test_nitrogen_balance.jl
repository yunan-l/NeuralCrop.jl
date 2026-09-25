using NeuralCrop
using Test

@testset "Daily nitrogen-balance diagnostics" begin
    crop = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    state = test_model_state(crop, soil)
    balance = @inferred init_nitrogen_balance(4, 1, identity)

    crop.state.nitrogen.total .= 1.0f0
    soil.nitrogen.nitrate .= 0.2f0
    soil.nitrogen.ammonium .= 0.1f0
    soil.nitrogen.litter .= 0.3f0
    soil.nitrogen.fast .= 0.4f0
    soil.nitrogen.slow .= 0.5f0
    NeuralCrop.record_nitrogen_balance_start!(balance, 1, state, state)
    NeuralCrop.record_nitrogen_balance_end!(balance, 1, state, state)

    @test balance.plant_before[1, 1] == 1.0f0
    @test balance.mineral_before[1, 1] ≈ 1.5f0 atol = 1.0f-6
    @test balance.organic_before[1, 1] ≈ 5.4f0 atol = 1.0f-6
    @test balance.residual[1, 1] ≈ 0.0f0 atol = 1.0f-6
    @test balance.relative_residual[1, 1] ≈ 0.0f0 atol = 1.0f-7

    # Root uptake is an internal soil-to-plant transfer, not a boundary input.
    NeuralCrop.record_nitrogen_balance_start!(balance, 2, state, state)
    soil.nitrogen.nitrate[1, 1] -= 0.25f0
    crop.state.nitrogen.total[1] += 0.25f0
    crop.fluxes.nitrogen.uptake[1] = 0.25f0
    NeuralCrop.record_nitrogen_balance_end!(balance, 2, state, state)

    @test balance.root_uptake[2, 1] == 0.25f0
    @test balance.residual[2, 1] ≈ 0.0f0 atol = 1.0f-6

    # Explicit external inputs must explain the corresponding stock increase.
    NeuralCrop.record_nitrogen_balance_start!(balance, 3, state, state)
    crop.state.nitrogen.total[1] += 0.7f0
    crop.fluxes.nitrogen.seed_input[1] = 0.7f0
    soil.nitrogen.ammonium[1, 1] += 0.2f0
    crop.fluxes.nitrogen.prescribed_fertilizer_input[1] = 0.2f0
    crop.fluxes.nitrogen.auto_fertilizer[1] = 0.3f0
    crop.state.nitrogen.total[1] += 0.3f0
    NeuralCrop.record_nitrogen_balance_end!(balance, 3, state, state)

    @test balance.seed_input[3, 1] == 0.7f0
    @test balance.prescribed_fertilizer_input[3, 1] == 0.2f0
    @test balance.automatic_fertilizer_input[3, 1] == 0.3f0
    @test balance.residual[3, 1] ≈ 0.0f0 atol = 2.0f-6

    # Harvest transfers residues internally and exports the remaining organ N.
    crop.state.nitrogen.total .= 1.0f0
    crop.state.nitrogen.leaf .= 0.2f0
    crop.state.nitrogen.root .= 0.3f0
    crop.state.nitrogen.storage .= 0.3f0
    crop.state.nitrogen.pool .= 0.2f0
    crop.fluxes.nitrogen.seed_input .= 0.0f0
    crop.fluxes.nitrogen.prescribed_fertilizer_input .= 0.0f0
    crop.fluxes.nitrogen.prescribed_manure_input .= 0.0f0
    crop.fluxes.nitrogen.auto_fertilizer .= 0.0f0
    crop.state.phenology.harvesting_previous .= false
    crop.state.phenology.harvesting .= true
    soil.nitrogen.litter .= 0.0f0
    soil.nitrogen.fast .= 0.0f0
    soil.nitrogen.slow .= 0.0f0
    soil.nitrogen.nitrate .= 0.0f0
    soil.nitrogen.ammonium .= 0.0f0
    soil.management.tillage_fraction .= Float32[1 0 0; 0 1 0; 0 0 1]
    output = init_output(1, identity)

    NeuralCrop.record_nitrogen_balance_start!(balance, 4, state, state)
    harvest_crop!(state, state, output, Float32[0.5], 100)
    # Follow the production daily loop: harvest_crop! sets is_growing = 0,
    # then crop_nitrogen! clears total and organ N through its inactive branch.
    crop_nitrogen!(
        state, cft1, state, zeros(Float32, 1), fill(20.0f0, 1);
        auto_fertilizer = false,
    )
    soil_nitrogen!(state, state)
    NeuralCrop.record_nitrogen_balance_end!(balance, 4, state, state)

    @test crop.state.nitrogen.total[1] == 0.0f0
    @test balance.harvest_export[4, 1] ≈ 0.5f0 atol = 1.0f-6
    @test balance.organic_after[4, 1] ≈ 0.5f0 atol = 1.0f-6
    @test balance.residual[4, 1] ≈ 0.0f0 atol = 1.0f-6
end

@testset "Atmospheric deposition nitrogen-balance input" begin
    crop = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    state = test_model_state(crop, soil)
    balance = init_nitrogen_balance(1, 1, identity)

    NeuralCrop.record_nitrogen_balance_start!(balance, 1, state, state)
    soil.nitrogen.nitrate[1, 1] += 0.11f0
    soil.nitrogen.ammonium[1, 1] += 0.09f0
    NeuralCrop.record_nitrogen_balance_end!(
        balance, 1, state, state;
        nitrate_deposition = Float32[0.11],
        ammonium_deposition = Float32[0.09],
    )

    @test balance.atmospheric_deposition[1, 1] == 0.2f0
    @test balance.residual[1, 1] ≈ 0.0f0 atol = 2.0f-6
end

@testset "Conservative baseline soil-N mineralization" begin
    crop = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    state = test_model_state(crop, soil)
    balance = init_nitrogen_balance(1, 1, identity)

    soil.nitrogen.litter .= 2.0f0
    soil.nitrogen.fast .= 3.0f0
    soil.nitrogen.slow .= 4.0f0
    soil.nitrogen.litter_response .= 0.1f0
    soil.decomposition.litter_response .= 1.0f0
    soil.decomposition.response .= 0.2f0
    # `init_soil` allocates state only. In a simulation `init_states!` loads
    # shared, normalized fast and slow c_shift distributions.
    soil.decomposition.shift_fast[1, 1] = 1.0f0
    soil.decomposition.shift_slow[1, 1] = 1.0f0

    NeuralCrop.record_nitrogen_balance_start!(balance, 1, state, state)
    soil_nitrogen!(state, state)
    NeuralCrop.record_nitrogen_balance_end!(balance, 1, state, state)

    @test balance.organic_after[1, 1] < balance.organic_before[1, 1]
    @test balance.mineral_after[1, 1] > balance.mineral_before[1, 1]
    @test balance.residual[1, 1] ≈ 0.0f0 atol = 1.0f-5
end
