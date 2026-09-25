using NeuralCrop
using Test

function warmup_array_snapshot(value)
    value isa AbstractArray && return Any[Array(value)]
    arrays = Any[]
    for name in fieldnames(typeof(value))
        append!(arrays, warmup_array_snapshot(getfield(value, name)))
    end
    return arrays
end

@testset "Agricultural warm-up" begin
    initial, _ = simulation_api_fixture(Float32)
    climate = climate_block(Float32, 365, 15, 1)
    simulation = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 3,
        diagnostics = true, fertilizer = :no,
    )
    initial_carbon = sum(Array(NeuralCrop.soil_carbon_prognostic(simulation.state).litter))

    report = agricultural_warmup!(simulation, climate)

    @test report.years == 10
    @test report.days == 3650
    @test report.forcing_years == 1
    @test length(report.consecutive_stable_years) == 1
    @test report.unconverged_cells == count(<(report.consecutive_years), report.consecutive_stable_years)
    @test size(report.soil.total_carbon) == (10, 1)
    @test size(report.soil.total_nitrogen) == (10, 1)
    @test size(report.soil.fast_carbon) == (10, 1)
    @test size(report.soil.slow_nitrogen) == (10, 1)
    @test all(isfinite, report.soil.total_carbon)
    @test all(report.soil.total_carbon .>= 0)
    @test all(report.soil.total_nitrogen .>= 0)
    @test all(report.soil.mineral_nitrogen .>= 0)
    @test size(report.calibrated_c_shift.fast) == (5, 1)
    @test size(report.calibrated_c_shift.slow) == (5, 1)
    @test all(report.calibrated_c_shift.response_sum .>= 0)
    @test sum(report.calibrated_c_shift.fast; dims = 1)[1] ≈ 1.0f0
    @test sum(report.calibrated_c_shift.slow; dims = 1)[1] ≈ 1.0f0
    @test all((0 .<= report.calibrated_pool_allocation.fast_carbon_fraction) .&
              (report.calibrated_pool_allocation.fast_carbon_fraction .<= 1))
    @test all((0 .<= report.calibrated_pool_allocation.fast_nitrogen_fraction) .&
              (report.calibrated_pool_allocation.fast_nitrogen_fraction .<= 1))
    drift = agricultural_warmup_drift(report)
    @test length(drift.total_carbon) == 11
    @test length(drift.total_nitrogen) == 11
    @test drift.spatial.cell_count == 1
    @test 0 <= drift.spatial.review_cell_fraction <= 1
    @test drift.recommendation in (:retain_40_60, :review_pool_allocation)
    @test report.soil.litter_carbon[end] != initial_carbon
    @test simulation.simulated_days == 0
    @test size(simulation.output.crop.npp, 1) == 0
    @test all(iszero, Array(simulation.carbon_balance.residual))

    constrained = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 3,
        diagnostics = false, fertilizer = :no,
    )
    constrained_carbon = NeuralCrop.soil_carbon_prognostic(constrained.state)
    constrained_nitrogen = NeuralCrop.soil_nitrogen_prognostic(constrained.state)
    target_carbon = Array(constrained_carbon.fast + constrained_carbon.slow)
    target_nitrogen = Array(
        constrained_nitrogen.fast + constrained_nitrogen.slow +
        constrained_nitrogen.nitrate + constrained_nitrogen.ammonium,
    )
    constrained_report = agricultural_warmup!(
        constrained, [climate];
        years = 2,
        maximum_years = 4,
        target_constrained = true,
        consecutive_years = 1,
        relative_tolerance = 1.0e6,
        pool_fraction_tolerance = 1.0e6,
        management_blocks = [(
            sdate = Int32[1],
            phu = Float32[-700],
            manure = Float32[0],
            fertilizer = Float32[0],
            residuefrac = Float32[0.5],
        )],
    )
    @test constrained_report.years == 2
    @test constrained_report.converged
    @test constrained_report.target_constrained
    @test constrained.state.inputs.crop.phenology.phu == Float32[700]
    @test constrained.state.inputs.crop.phenology.winter_type == Bool[true]
    @test constrained_report.converged_cell_fraction == 1
    @test size(constrained_report.target_correction.carbon) == (2, 1)
    constrained_drift = agricultural_warmup_drift(constrained_report)
    @test constrained_drift.recommendation == :target_constrained_converged
    @test constrained_drift.convergence.actual_years == 2
    @test Array(constrained_carbon.fast + constrained_carbon.slow) ≈ target_carbon
    @test Array(
        constrained_nitrogen.fast + constrained_nitrogen.slow +
        constrained_nitrogen.nitrate + constrained_nitrogen.ammonium,
    ) ≈ target_nitrogen

    stable_snapshot = (
        total_carbon = Float32[100], total_nitrogen = Float32[10],
        fast_carbon = Float32[40], slow_carbon = Float32[60],
        fast_nitrogen = Float32[4], slow_nitrogen = Float32[6],
    )
    stable_correction = (carbon = Float32[-5], nitrogen = Float32[-2])
    consecutive = [0]
    @test NeuralCrop._warmup_convergence!(
        consecutive, stable_snapshot, stable_snapshot, stable_snapshot,
        stable_correction, stable_correction;
        relative_tolerance = 0.01,
        pool_fraction_tolerance = 0.01,
        consecutive_years = 1,
    ) == 1
    @test consecutive == [1]

    reduced_consecutive = [0]
    @test NeuralCrop._warmup_convergence!(
        reduced_consecutive, stable_snapshot, stable_snapshot, stable_snapshot,
        stable_correction, stable_correction;
        relative_tolerance = 0.01,
        pool_fraction_tolerance = 0.01,
        consecutive_years = 1,
        convergence_reducer = (converged, total) -> (converged + 2, total + 3),
    ) == 0.75
    @test reduced_consecutive == [1]

    reduced = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 3,
        diagnostics = false, fertilizer = :no,
    )
    reduced_report = agricultural_warmup!(
        reduced, [climate];
        years = 2,
        maximum_years = 2,
        target_constrained = true,
        consecutive_years = 1,
        relative_tolerance = 1.0e6,
        pool_fraction_tolerance = 1.0e6,
        required_converged_fraction = 0.75,
        convergence_reducer = (converged, total) -> (converged + 2, total + 3),
    )
    @test reduced_report.converged
    @test reduced_report.converged_cell_fraction == 0.75
    @test reduced_report.unconverged_cells == 1

    capped = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 3,
        diagnostics = false, fertilizer = :no,
    )
    capped_report = agricultural_warmup!(
        capped, [climate];
        years = 2,
        maximum_years = 3,
        target_constrained = true,
        consecutive_years = 4,
    )
    @test capped_report.years == 3
    @test !capped_report.converged
    @test capped_report.unconverged_cells == 1
    @test agricultural_warmup_drift(capped_report).recommendation ==
        :target_constrained_maximum_years

    two_year_climate = climate_block(Float32, 730, 15, 1)
    incomplete_cycle = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 3,
        diagnostics = false, fertilizer = :no,
    )
    incomplete_report = agricultural_warmup!(
        incomplete_cycle, two_year_climate;
        years = 1,
        maximum_years = 1,
        target_constrained = true,
        consecutive_years = 1,
        relative_tolerance = 1.0e6,
        pool_fraction_tolerance = 1.0e6,
    )
    @test incomplete_report.forcing_years == 2
    @test !incomplete_report.converged
    @test incomplete_report.converged_cell_fraction == 0

    complete_cycle = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 3,
        diagnostics = false, fertilizer = :no,
    )
    complete_report = agricultural_warmup!(
        complete_cycle, two_year_climate;
        years = 1,
        maximum_years = 3,
        target_constrained = true,
        consecutive_years = 1,
        relative_tolerance = 1.0e6,
        pool_fraction_tolerance = 1.0e6,
    )
    @test complete_report.years == 3
    @test complete_report.converged

    run_simulation!(simulation, climate; end_day = 3, spinup = false)
    @test simulation.simulated_days == 3
    @test size(simulation.output.crop.npp) == (3, 1)

    started = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 3,
        diagnostics = false, fertilizer = :no,
    )
    transition_day!(started, climate)
    @test_throws ArgumentError agricultural_warmup!(started, climate)
    @test_throws ArgumentError agricultural_warmup!(
        initialize_simulation(
            cft1, initial;
            indices = [1], T = Float32, days = 3,
            diagnostics = false, fertilizer = :no,
        ),
        climate_block(Float32, 364, 15, 1),
    )

    eager_climate = (
        temp = reshape(Float32[15 + 2 * (day > 365) for day in 1:730], 730, 1),
        prec = reshape(Float32[1 + (day % 11 == 0) for day in 1:730], 730, 1),
        sw = fill(180.0f0, 730, 1),
        lw = fill(-40.0f0, 730, 1),
        co2 = reshape(Float32[400 + (day > 365) for day in 1:730], 730, 1),
        co2_daily = true,
        backend_neutral = true,
    )
    ranges = [1:73, 74:289, 290:410, 411:619, 620:730]
    streamed_climate = [(
        temp = eager_climate.temp[range, :],
        prec = eager_climate.prec[range, :],
        sw = eager_climate.sw[range, :],
        lw = eager_climate.lw[range, :],
        co2 = eager_climate.co2[range, :],
        co2_daily = true,
        backend_neutral = true,
    ) for range in ranges]
    eager_simulation = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 3,
        diagnostics = true, fertilizer = :no,
    )
    streamed_simulation = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 3,
        diagnostics = true, fertilizer = :no,
    )
    eager_report = agricultural_warmup!(eager_simulation, eager_climate; years = 3)
    streamed_report = agricultural_warmup!(streamed_simulation, streamed_climate; years = 3)

    @test streamed_report == eager_report
    @test warmup_array_snapshot(streamed_simulation.state.prognostic) ==
        warmup_array_snapshot(eager_simulation.state.prognostic)
    @test streamed_simulation.simulated_days == 0
    @test size(streamed_simulation.output.crop.npp, 1) == 0
    @test all(iszero, Array(streamed_simulation.carbon_balance.residual))

    mktempdir() do directory
        checkpoint_path = joinpath(directory, "warmup_checkpoint.jld2")
        @test save_checkpoint(checkpoint_path, streamed_simulation) == checkpoint_path
        restored = initialize_simulation(
            cft1, initial;
            indices = [1], T = Float32, days = 3,
            diagnostics = true, fertilizer = :no,
        )
        @test restore_checkpoint!(restored, checkpoint_path) === restored
        @test restored.simulated_days == 0
        @test warmup_array_snapshot(restored.state.prognostic) ==
            warmup_array_snapshot(streamed_simulation.state.prognostic)

        run_simulation!(streamed_simulation, streamed_climate[1]; end_day = 3, spinup = false)
        run_simulation!(restored, streamed_climate[1]; end_day = 3, spinup = false)
        @test restored.output.crop.npp == streamed_simulation.output.crop.npp
        @test warmup_array_snapshot(restored.state.prognostic) ==
            warmup_array_snapshot(streamed_simulation.state.prognostic)
    end

    @test_throws ArgumentError agricultural_warmup!(
        initialize_simulation(
            cft1, initial;
            indices = [1], T = Float32, days = 3,
            diagnostics = false, fertilizer = :no,
        ),
        streamed_climate[1:2],
    )
end
