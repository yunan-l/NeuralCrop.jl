using NeuralCrop
using Test

@testset "Surface-litter water capacity and interception" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    state = test_model_state(soil)
    litter_carbon = 20.0f0 * 0.42f0 * 71.1f0
    soil.carbon.litter[1, 1] = litter_carbon

    update_surface_litter_properties!(state)

    expected_dry_matter = litter_carbon / 0.42f0
    expected_capacity = 2.0f-3 * expected_dry_matter
    expected_cover = 1.0f0 - exp(-6.0f-3 * expected_dry_matter)
    @test soil.surface_litter.depth[1] ≈ 0.02f0 atol = 1.0f-6
    @test soil.surface_litter.water_capacity[1] ≈ expected_capacity
    @test soil.surface_litter.cover[1] ≈ expected_cover

    soil.water.infiltration .= 10.0f0
    surface_litter_interception!(state)
    @test soil.surface_litter.interception[1] ≈ expected_capacity
    @test soil.surface_litter.water_storage[1] ≈ expected_capacity
    @test soil.water.infiltration[1] ≈ 10.0f0 - expected_capacity
end

@testset "Reduced litter capacity conserves water" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    state = test_model_state(soil)
    soil.carbon.litter[1, 1] = 100.0f0
    update_surface_litter_properties!(state)
    soil.surface_litter.water_storage .= soil.surface_litter.water_capacity
    soil.water.storage[1, 1] = 5.0f0
    water_before = soil.water.storage[1, 1] +
                   soil.surface_litter.water_storage[1]

    soil.carbon.litter[1, 1] = 50.0f0
    update_surface_litter_properties!(state)

    water_after = soil.water.storage[1, 1] +
                  soil.surface_litter.water_storage[1]
    @test water_after ≈ water_before atol = 1.0f-6
    @test soil.surface_litter.water_storage[1] ≈
          soil.surface_litter.water_capacity[1]
end

@testset "Wet litter evaporation" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)
    state = test_model_state(crop, soil)
    soil.carbon.litter[1, 1] = 20.0f0 * 0.42f0 * 71.1f0
    update_surface_litter_properties!(state)
    soil.surface_litter.water_storage .= soil.surface_litter.water_capacity
    storage_before = soil.surface_litter.water_storage[1]

    soil.water.storage .= 50.0f0
    soil.water.wilting_storage .= 10.0f0
    soil.water.holding_capacity_storage .= 100.0f0
    evaporation!(Float32[2.0], state, state)

    @test soil.surface_litter.evaporation[1] > 0.0f0
    @test soil.surface_litter.water_storage[1] ≈
          storage_before - soil.surface_litter.evaporation[1]
    @test soil.surface_litter.water_storage[1] >= 0.0f0
end
