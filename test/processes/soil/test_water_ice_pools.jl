using NeuralCrop
using Test

@testset "LPJmL three-reservoir water-ice partition" begin
    # 20 mm below PWP + 40 mm plant-available water; half is frozen.
    pwp_ice, available_ice, free_ice, relative_water, free_water =
        NeuralCrop.lpjml_water_ice_partition(60.0f0, 30.0f0, 20.0f0, 40.0f0)
    @test pwp_ice ≈ 0.5f0
    @test available_ice ≈ 20.0f0
    @test free_ice == 0.0f0
    @test relative_water ≈ 0.5f0
    @test free_water == 0.0f0

    # Free water remains liquid until the PWP and available pools are frozen.
    pwp_ice, available_ice, free_ice, relative_water, free_water =
        NeuralCrop.lpjml_water_ice_partition(80.0f0, 40.0f0, 20.0f0, 40.0f0)
    @test pwp_ice ≈ 2.0f0 / 3.0f0
    @test available_ice ≈ 80.0f0 / 3.0f0
    @test free_ice == 0.0f0
    @test relative_water ≈ 1.0f0 / 3.0f0
    @test free_water == 20.0f0

    # Free ice appears only after both bound reservoirs are completely frozen.
    pwp_ice, available_ice, free_ice, relative_water, free_water =
        NeuralCrop.lpjml_water_ice_partition(80.0f0, 70.0f0, 20.0f0, 40.0f0)
    @test pwp_ice == 1.0f0
    @test available_ice == 40.0f0
    @test free_ice == 10.0f0
    @test relative_water == 0.0f0
    @test free_water == 10.0f0
end

@testset "Layer states retain LPJmL ice-pool invariants" begin
    soil = init_soil(2, soilparams.soildepth, identity)
    state = test_model_state(soil)
    soil.water.storage .= 50.0f0
    pedotransfer!(state)
    soil_temperature!(state, Float32[-20.0, 10.0], Float32[2.0, 2.0])

    component_ice =
        soil.water.wilting_ice_fraction .* soil.water.wilting_storage .+
        soil.water.available_ice_storage .+
        soil.water.free_ice_storage
    @test component_ice ≈ soil.water.ice_storage atol = 2.0f-5
    bound_not_fully_frozen = soil.water.wilting_ice_fraction .< 1.0f0
    @test all(soil.water.free_ice_storage[bound_not_fully_frozen] .== 0.0f0)
    @test all((0.0f0 .<= soil.water.relative_content) .&
              (soil.water.relative_content .<= 1.0f0))
end
