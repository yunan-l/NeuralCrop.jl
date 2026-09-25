using NeuralCrop
using Test

@testset "Staged daily soil-water update" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)
    state = test_model_state(crop, soil)

    soil.properties.sand_fraction .= 0.4f0
    soil.properties.clay_fraction .= 0.2f0
    soil.water.storage .= Float32[40, 60, 100, 200, 200]
    soil.nitrogen.nitrate .= reshape(Float32[10, 20, 30, 40, 50], 5, 1)
    crop.fluxes.water.interception .= 0.0f0
    crop.auxiliary.root.distribution .= 0.0f0
    crop.auxiliary.root.distribution[1] = 1.0f0
    pedotransfer!(state)

    storage_before = sum(soil.water.storage)
    nitrate_before = sum(soil.nitrogen.nitrate)
    water_availability_before = sum(soil.water.relative_content .* crop.auxiliary.root.distribution)
    soil_infiltration!(state, state, Float32[10.0])
    storage_after_infiltration = sum(soil.water.storage)
    water_availability_after = sum(soil.water.relative_content .* crop.auxiliary.root.distribution)

    @test storage_after_infiltration > storage_before
    @test storage_after_infiltration + soil.water.surface_runoff[1] +
          sum(soil.water.lateral_runoff) + soil.water.bottom_drainage[1] ≈ storage_before + 10.0f0 atol = 1.0f-4
    @test water_availability_after >= water_availability_before
    @test sum(soil.nitrogen.nitrate) + soil.nitrogen.leaching[1] ≈ nitrate_before atol = 1.0f-4

    crop.fluxes.water.transpiration_layer .= 0.4f0
    soil.water.evaporation .= 0.2f0
    soil_evapotranspiration!(state, state)

    @test sum(soil.water.storage) ≈ storage_after_infiltration - 3.0f0 atol = 1.0f-5
end

@testset "Full irrigation restores total water to field capacity" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)
    state = test_model_state(crop, soil)

    soil.properties.sand_fraction .= 0.4f0
    soil.properties.clay_fraction .= 0.2f0
    pedotransfer!(state)

    target_storage = soil.water.field_capacity .* reshape(soil.properties.layer_depth, :, 1)

    soil.water.storage .= 0.0f0
    soil.water.ice_storage .= 0.0f0
    soil_evapotranspiration!(state, state; irrigation = true)
    @test soil.water.storage == target_storage

    ice_storage = min.(target_storage .* 0.25f0, 10.0f0)
    soil.water.storage .= 0.0f0
    soil.water.ice_storage .= ice_storage
    soil_evapotranspiration!(state, state; irrigation = true)

    @test soil.water.storage .+ soil.water.ice_storage ≈ target_storage
    @test soil.water.storage ≈ target_storage .- ice_storage
end

@testset "LPJmL infiltration iteration cap preserves water" begin
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)
    state = test_model_state(crop, soil)
    soil.properties.sand_fraction .= 0.4f0
    soil.properties.clay_fraction .= 0.2f0
    soil.water.storage .= Float32[40, 60, 100, 200, 200]
    crop.fluxes.water.interception .= 0.0f0
    pedotransfer!(state)

    storage_before = sum(soil.water.storage)
    # 1001 four-millimetre slugs force LPJmL's MAXITER fallback.
    precipitation = 4004.0f0
    soil_infiltration!(state, state, Float32[precipitation])
    accounted_water = sum(soil.water.storage) + soil.water.surface_runoff[1] +
        sum(soil.water.lateral_runoff) + soil.water.bottom_drainage[1]
    @test accounted_water ≈ storage_before + precipitation atol = 0.2f0
end

@testset "Same-day rain pulse affects crop water stress" begin
    dry_soil = init_soil(1, soilparams.soildepth, identity)
    dry_crop = init_crop(1, identity)
    pet = init_pet(1, identity)
    dry_state = test_model_state(dry_crop, dry_soil; pet)

    dry_soil.properties.sand_fraction .= 0.4f0
    dry_soil.properties.clay_fraction .= 0.2f0
    dry_soil.water.storage .= Float32[25, 60, 100, 200, 200]
    dry_crop.state.phenology.is_growing .= 1
    dry_crop.state.carbon.root .= 50.0f0
    dry_crop.auxiliary.root.distribution .= 0.0f0
    dry_crop.auxiliary.root.distribution[1] = 1.0f0
    dry_crop.auxiliary.canopy.canopy_wet .= 0.0f0
    dry_crop.fluxes.water.interception .= 0.0f0
    pet.eeq .= 5.0f0
    pet.daylength .= 12.0f0
    pedotransfer!(dry_state)

    wet_soil = deepcopy(dry_soil)
    wet_crop = deepcopy(dry_crop)
    wet_state = test_model_state(wet_crop, wet_soil; pet)
    soil_infiltration!(wet_state, wet_state, Float32[20.0])

    assimilation = Float32[5.0]
    co2 = Float32[40.0]
    transpiration!(assimilation, cft1, dry_state, pet, dry_state, co2)
    transpiration!(assimilation, cft1, wet_state, pet, wet_state, co2)

    expected_rootzone_water = sum(
        wet_soil.water.relative_content[1:3, 1] .*
        wet_soil.water.holding_capacity_storage[1:3, 1] .*
        wet_crop.auxiliary.root.distribution[1:3],
    )
    @test wet_crop.state.water.supply_sum[1] > dry_crop.state.water.supply_sum[1]
    @test wet_crop.state.water.sufficiency[1] > dry_crop.state.water.sufficiency[1]
    @test sum(wet_crop.fluxes.water.transpiration_layer) > sum(dry_crop.fluxes.water.transpiration_layer)
    @test wet_crop.auxiliary.root.zone_available_water[1] ≈ expected_rootzone_water atol = 1.0f-6
end
