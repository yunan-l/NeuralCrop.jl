using NeuralCrop
using Test

@testset "Daily freeze-thaw diagnostics" begin
    soil = init_soil(2, soilparams.soildepth, identity)
    state = test_model_state(soil)
    diagnostics = init_thermal_balance(1, 2, identity)
    soil.water.storage .= 50.0f0
    pedotransfer!(state)
    soil_temperature!(state, Float32[-20.0, 10.0], Float32[2.0, 2.0])
    NeuralCrop.record_thermal_balance!(diagnostics, 1, state)

    @test diagnostics.total_ice_storage[1, 1] > 0.0f0
    @test diagnostics.total_ice_storage[1, 2] == 0.0f0
    @test diagnostics.maximum_frozen_fraction[1, 1] > 0.0f0
    @test diagnostics.wilting_ice_storage[1, 1] > 0.0f0
    @test diagnostics.available_ice_storage[1, 1] > 0.0f0
    @test diagnostics.free_ice_storage[1, 1] >= 0.0f0
    @test diagnostics.ice_pool_residual[1, 1] <= 2.0f-5
    @test all(isfinite, diagnostics.energy_residual)
    @test all(iszero, diagnostics.untracked_water_energy_flux)
    @test all(isfinite, diagnostics.percolation_energy_residual)
    @test all(iszero, diagnostics.rain_energy_input)
    @test all(iszero, diagnostics.snowmelt_energy_input)
    @test all(iszero, diagnostics.lateral_runoff_energy_output)
    @test all(iszero, diagnostics.bottom_drainage_energy_output)
    @test diagnostics.minimum_temperature[1, 1] <=
          diagnostics.maximum_temperature[1, 1]
end
