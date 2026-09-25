using NeuralCrop
using Test

@testset "Five-layer implicit soil heat conduction" begin
    soil = init_soil(2, soilparams.soildepth, identity)
    state = test_model_state(soil)
    soil.thermal.diffusivity_0 .= 0.6f0
    soil.thermal.diffusivity_15 .= 0.7f0
    soil.water.relative_content .= 0.1f0

    # A profile in equilibrium with its upper boundary remains unchanged.
    soil_temperature!(state, Float32[10.0, 10.0], Float32[10.0, 10.0])
    @test all(soil.thermal.initialized)
    @test all(isapprox.(soil.thermal.temperature, 10.0f0; atol = 1.0f-5))

    # A surface warming step creates a finite, depth-damped vertical gradient.
    soil_temperature!(state, Float32[30.0, -10.0], Float32[10.0, 10.0])
    warm_profile = soil.thermal.temperature[:, 1]
    cold_profile = soil.thermal.temperature[:, 2]

    @test all(isfinite, soil.thermal.temperature)
    @test 10.0f0 < warm_profile[1] < 30.0f0
    @test all(diff(warm_profile) .< 0.0f0)
    @test -10.0f0 < cold_profile[1] < 10.0f0
    @test all(diff(cold_profile) .> 0.0f0)
    @test warm_profile[1] - 10.0f0 > warm_profile[5] - 10.0f0
    @test 10.0f0 - cold_profile[1] > 10.0f0 - cold_profile[5]

    # Backward Euler remains bounded under sustained forcing.
    for _ in 1:100
        soil_temperature!(state, Float32[30.0, -10.0])
    end
    @test all((soil.thermal.temperature[:, 1] .>= 10.0f0) .&
              (soil.thermal.temperature[:, 1] .<= 30.0f0))
    @test all((soil.thermal.temperature[:, 2] .>= -10.0f0) .&
              (soil.thermal.temperature[:, 2] .<= 10.0f0))
end

@testset "Snow thermal resistance" begin
    bare_soil = init_soil(1, soilparams.soildepth, identity)
    snow_soil = init_soil(1, soilparams.soildepth, identity)
    bare_state = test_model_state(bare_soil)
    snow_state = test_model_state(snow_soil)
    for soil in (bare_soil, snow_soil)
        soil.thermal.diffusivity_0 .= 0.6f0
        soil.thermal.diffusivity_15 .= 0.7f0
        soil.water.relative_content .= 0.1f0
    end
    snow_soil.snow.height .= 0.67f0

    soil_temperature!(bare_state, Float32[30.0], Float32[0.0])
    soil_temperature!(snow_state, Float32[30.0], Float32[0.0])

    @test 0.0f0 < snow_soil.thermal.temperature[1, 1] <
                  bare_soil.thermal.temperature[1, 1]
    @test sum(snow_soil.thermal.temperature) <
          sum(bare_soil.thermal.temperature)
    @test all(isfinite, snow_soil.thermal.temperature)

    # Snow also damps a cold surface pulse rather than changing its direction.
    bare_soil.thermal.temperature .= 0.0f0
    snow_soil.thermal.temperature .= 0.0f0
    soil_temperature!(bare_state, Float32[-30.0])
    soil_temperature!(snow_state, Float32[-30.0])
    @test bare_soil.thermal.temperature[1, 1] <
          snow_soil.thermal.temperature[1, 1] < 0.0f0
end

@testset "Dry surface-litter thermal resistance" begin
    bare_soil = init_soil(1, soilparams.soildepth, identity)
    litter_soil = init_soil(1, soilparams.soildepth, identity)
    bare_state = test_model_state(bare_soil)
    litter_state = test_model_state(litter_soil)
    for soil in (bare_soil, litter_soil)
        soil.thermal.diffusivity_0 .= 0.6f0
        soil.thermal.diffusivity_15 .= 0.7f0
        soil.water.relative_content .= 0.1f0
    end

    # LPJmL test convention: this carbon stock corresponds to 2 cm dry litter.
    litter_soil.carbon.litter[1, 1] =
        20.0f0 * soil_thermal_params.litter_carbon_fraction *
        soil_thermal_params.litter_bulk_density
    update_surface_litter_properties!(litter_state)
    @test litter_soil.surface_litter.depth[1] ≈ 0.02f0 atol = 1.0f-6

    soil_temperature!(bare_state, Float32[30.0], Float32[0.0])
    soil_temperature!(litter_state, Float32[30.0], Float32[0.0])
    @test 0.0f0 < litter_soil.thermal.temperature[1, 1] <
                  bare_soil.thermal.temperature[1, 1]
    @test sum(litter_soil.thermal.temperature) <
          sum(bare_soil.thermal.temperature)

    bare_soil.thermal.temperature .= 0.0f0
    litter_soil.thermal.temperature .= 0.0f0
    soil_temperature!(bare_state, Float32[-30.0])
    soil_temperature!(litter_state, Float32[-30.0])
    @test bare_soil.thermal.temperature[1, 1] <
          litter_soil.thermal.temperature[1, 1] < 0.0f0
end

@testset "Wet surface litter conducts more heat than dry litter" begin
    dry_soil = init_soil(1, soilparams.soildepth, identity)
    wet_soil = init_soil(1, soilparams.soildepth, identity)
    dry_state = test_model_state(dry_soil)
    wet_state = test_model_state(wet_soil)
    for soil in (dry_soil, wet_soil)
        soil.thermal.diffusivity_0 .= 0.6f0
        soil.thermal.diffusivity_15 .= 0.7f0
        soil.water.relative_content .= 0.1f0
        soil.carbon.litter[1, 1] =
            20.0f0 * soil_thermal_params.litter_carbon_fraction *
            soil_thermal_params.litter_bulk_density
    end
    update_surface_litter_properties!(dry_state)
    update_surface_litter_properties!(wet_state)
    wet_soil.surface_litter.water_storage .=
        wet_soil.surface_litter.water_capacity

    soil_temperature!(dry_state, Float32[30.0], Float32[0.0])
    soil_temperature!(wet_state, Float32[30.0], Float32[0.0])
    @test wet_soil.surface_litter.conductivity[1] >
          dry_soil.surface_litter.conductivity[1]
    @test dry_soil.thermal.temperature[1, 1] <
          wet_soil.thermal.temperature[1, 1]
    @test isfinite(wet_soil.surface_litter.temperature[1])
end

@testset "LPJmL soil decomposition temperature bounds" begin
    temperatures = Float32[-30.0, -15.01, -15.0, 10.0, 40.0, 55.0]
    response = NeuralCrop.temp_response(temperatures)

    @test response[1] == 0.0f0
    @test response[2] == 0.0f0
    @test response[3] > 0.0f0
    @test response[4] ≈ 1.0f0 atol = 1.0f-6
    @test response[5] == response[6]
    @test all(isfinite, response)
end
