using NeuralCrop
using Test

function initialized_thermal_soil(; storage = Float32[40, 60, 100, 200, 200],
                                  temperature = 5.0f0)
    soil = init_soil(1, soilparams.soildepth, identity)
    crop = init_crop(1, identity)
    state = test_model_state(crop, soil)
    soil.properties.sand_fraction .= 0.4f0
    soil.properties.clay_fraction .= 0.2f0
    soil.water.storage .= reshape(storage, 5, 1)
    crop.fluxes.water.interception .= 0.0f0
    pedotransfer!(state)
    soil_temperature!(state, Float32[temperature], Float32[temperature])
    return soil, crop, state
end

@testset "LPJmL percolation enthalpy" begin
    @testset "dry upper boundary has no rain or melt energy" begin
        soil, crop, state = initialized_thermal_soil()
        soil_infiltration!(
            state, state, Float32[0.0];
            snowmelt = Float32[0.0], air_temperature = Float32[5.0],
        )

        @test iszero(soil.thermal.rain_energy_input[1])
        @test iszero(soil.thermal.snowmelt_energy_input[1])
        @test all(isfinite, soil.thermal.temperature)
        @test all(isfinite, soil.thermal.enthalpy)
    end

    @testset "rain at soil temperature preserves temperature" begin
        soil, crop, state = initialized_thermal_soil()
        soil_infiltration!(
            state, state, Float32[2.0];
            snowmelt = Float32[0.0], air_temperature = Float32[5.0],
        )

        accepted_rain = soil.water.influx[1, 1]
        expected_energy = accepted_rain * 0.001f0 *
            (soil_thermal_params.volumetric_fusion_heat +
             soil_thermal_params.water_heat_capacity * 5.0f0)
        @test soil.thermal.rain_energy_input[1] ≈ expected_energy rtol = 2.0f-6
        @test soil.thermal.snowmelt_energy_input[1] == 0.0f0
        @test maximum(abs.(soil.thermal.temperature .- 5.0f0)) < 2.0f-5
        @test abs(soil.thermal.percolation_energy_residual[1]) < 1.0f0
        @test all(iszero, soil.thermal.percolation_energy)
    end

    @testset "rain and meltwater retain separate upper-boundary enthalpy" begin
        soil, crop, state = initialized_thermal_soil(; temperature = 10.0f0)
        soil_infiltration!(
            state, state, Float32[4.0];
            snowmelt = Float32[2.0], air_temperature = Float32[10.0],
        )

        accepted = soil.water.influx[1, 1]
        expected_rain = accepted * 0.5f0 * 0.001f0 *
            (soil_thermal_params.volumetric_fusion_heat +
             soil_thermal_params.water_heat_capacity * 10.0f0)
        expected_melt = accepted * 0.5f0 * 0.001f0 *
            soil_thermal_params.volumetric_fusion_heat
        @test soil.thermal.rain_energy_input[1] ≈ expected_rain rtol = 2.0f-6
        @test soil.thermal.snowmelt_energy_input[1] ≈ expected_melt rtol = 2.0f-6
        @test soil.thermal.temperature[1, 1] < 10.0f0
        @test abs(soil.thermal.percolation_energy_residual[1]) < 1.0f0
    end

    @testset "layer transfers cancel and boundary ledger closes" begin
        soil, crop, state = initialized_thermal_soil(
            storage = Float32[160, 240, 480, 800, 1600],
            temperature = 10.0f0,
        )
        soil_infiltration!(
            state, state, Float32[50.0];
            snowmelt = Float32[10.0], air_temperature = Float32[10.0],
        )

        @test any(soil.water.outflux .> 0.0f0)
        boundary_scale = max(
            abs(soil.thermal.rain_energy_input[1]) +
            abs(soil.thermal.snowmelt_energy_input[1]) +
            abs(soil.thermal.lateral_runoff_energy_output[1]) +
            abs(soil.thermal.bottom_drainage_energy_output[1]),
            1.0f0,
        )
        @test abs(soil.thermal.percolation_energy_residual[1]) / boundary_scale < 5.0f-6
        @test all(isfinite, soil.thermal.temperature)
        @test all(isfinite, soil.thermal.enthalpy)
    end
end
