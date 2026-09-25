using NeuralCrop
using JLD2
using Test

function simulation_api_fixture(::Type{T}) where {T <: AbstractFloat}
    cells = 1
    layers = 5
    initial = (
        coords = [1],
        latitude = T[45],
        crop = (
            sdate = Int32[1],
            phu = T[543],
            manure = zeros(T, cells),
            fertilizer = T[24.55],
            residuefrac = T[0.67],
        ),
        soilparam = (
            soilph = T[6.5],
            w_sat = fill(T(0.45), layers, cells),
            sand = T[0.4],
            clay = T[0.2],
            tdiff_0 = T[0.7],
            tdiff_15 = T[0.75],
            soildepth = T[200, 300, 500, 1000, 1000],
        ),
        initial_state = (
            swc = reshape(T[57.41, 55.32, 126.13, 274.59, 285.71], layers, cells),
            litc = reshape(T[0.13, 187.5, 225.36], 3, cells),
            fastc = fill(T(10), layers, cells),
            slowc = fill(T(100), layers, cells),
            litn = reshape(T[0.0047, 6.47, 9.47], 3, cells),
            fastn = fill(T(1), layers, cells),
            slown = fill(T(10), layers, cells),
        ),
    )
    climate = (
        temp_spinup = fill(T(10), 365, cells),
        temp = fill(T(15), 3, cells),
        prec = fill(T(1), 3, cells),
        swdown = fill(T(180), 3, cells),
        lwnet = fill(T(-40), 3, cells),
        windspeed = fill(T(2), 3, cells),
        co2 = T[400],
    )
    return initial, climate
end

function climate_block(::Type{T}, days, temperature, precipitation) where {T <: AbstractFloat}
    annual_co2 = fill(T(400), cld(days, 365))
    return (
        temp_spinup = fill(T(10), 365, 1),
        temp = fill(T(temperature), days, 1),
        prec = fill(T(precipitation), days, 1),
        swdown = fill(T(180), days, 1),
        lwnet = fill(T(-40), days, 1),
        windspeed = fill(T(2), days, 1),
        co2 = annual_co2,
    )
end

function runtime_array_ids(value)
    value isa AbstractArray && return UInt[objectid(value)]
    ids = UInt[]
    for name in fieldnames(typeof(value))
        append!(ids, runtime_array_ids(getfield(value, name)))
    end
    return ids
end

@testset "Simulation memory estimate" begin
    initial, _ = simulation_api_fixture(Float32)
    simulation = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 4,
        diagnostics = false, fertilizer = :no,
    )
    estimate = estimate_memory(simulation; block_days = 2, safety_factor = 1)
    preallocation = estimate_memory(
        1, 4; T = Float32, diagnostics = false, block_days = 2,
        backend = :cpu, safety_factor = 1,
    )
    prefetched = estimate_memory(
        simulation; block_days = 2, prefetch = true, safety_factor = 1,
    )

    @test estimate.cells == 1
    @test estimate.days == 4
    @test estimate.backend == :cpu
    @test estimate.diagnostics_bytes == 0
    @test estimate.projected_output_bytes == 356
    @test preallocation == estimate
    @test estimate_memory(
        1, 1; T = Float64, diagnostics = false, block_days = 1,
        backend = :cpu, safety_factor = 1,
    ).persistent_state_bytes == 13342
    @test estimate.forcing_block_bytes == 56
    @test prefetched.host_forcing_bytes == estimate.host_forcing_bytes + 56
    @test prefetched.host_peak_bytes == estimate.host_peak_bytes + 56
    warmup_estimate = estimate_memory(
        simulation; block_days = 2, warmup_years = 10, safety_factor = 1,
    )
    @test warmup_estimate.warmup_history_bytes == 10 * 12 * sizeof(Float32)
    @test warmup_estimate.host_peak_bytes ==
        estimate.host_peak_bytes + warmup_estimate.warmup_history_bytes
    cached_estimate = estimate_memory(
        simulation; block_days = 2, cached_forcing_blocks = 3, safety_factor = 1,
    )
    @test cached_estimate.cached_forcing_bytes == 3 * estimate.forcing_block_bytes
    @test cached_estimate.host_peak_bytes ==
        estimate.host_peak_bytes + cached_estimate.cached_forcing_bytes
    @test estimate.device_peak_bytes == 0
    @test estimate.recommended_host_peak_gib == estimate.host_peak_bytes / 2.0^30
    stream = OutputStream(
        [OutputVariable(:crop, :npp)]; cell_ids = Int32[1],
    )
    streamed = estimate_memory(
        simulation; block_days = 2, output_stream = stream, safety_factor = 1,
    )
    @test streamed.streaming_output
    @test streamed.projected_output_bytes < estimate.projected_output_bytes
    @test streamed.output_growth_bytes == 0
    reader = (block_days = 2, selection = (compact_indices = [1],))
    @test estimate_memory(simulation, reader; safety_factor = 1) == estimate
    mismatched = (block_days = 2, selection = (compact_indices = [1, 2],))
    @test_throws DimensionMismatch estimate_memory(simulation, mismatched)
    @test_throws ArgumentError estimate_memory(simulation; block_days = 0)
    @test_throws ArgumentError estimate_memory(0, 4; block_days = 2)
    @test_throws ArgumentError estimate_memory(1, 4; block_days = 2, backend = :gpu)
    @test_throws ArgumentError estimate_memory(
        simulation; block_days = 2, safety_factor = 0.9,
    )
    @test_throws ArgumentError estimate_memory(
        simulation; block_days = 2, warmup_years = -1,
    )

    with_diagnostics = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 4,
        diagnostics = true, fertilizer = :no,
    )
    allocated = estimate_memory(with_diagnostics; block_days = 2, safety_factor = 1)
    predicted = estimate_memory(
        1, 4; T = Float32, diagnostics = true, block_days = 2,
        backend = :cpu, safety_factor = 1,
    )
    @test predicted == allocated
end

@testset "Nitrogen limitation remains explicit" begin
    default_configuration = SimulationConfiguration(
        Float32, identity, 1, Int64[1], Int32[1],
    )
    @test !default_configuration.nitrogen_limit_vcmax

    configuration = SimulationConfiguration(
        Float32, identity, 1, Int64[1], Int32[1];
        nitrogen_limit_vcmax = true,
    )
    @test configuration.nitrogen_limit_vcmax
end

@testset "Crop respiration mode remains explicit" begin
    default_configuration = SimulationConfiguration(
        Float32, identity, 1, Int64[1], Int32[1],
    )
    @test !default_configuration.crop_resp_fix

    configuration = SimulationConfiguration(
        Float32, identity, 1, Int64[1], Int32[1];
        crop_resp_fix = true,
    )
    @test configuration.crop_resp_fix
end

@testset "Nitrogen-limited daily driver accepts deposition forcing" begin
    initial, climate = simulation_api_fixture(Float32)
    simulation = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 3, fertilizer = :no,
        nitrogen_limit_vcmax = true,
    )
    deposition_climate = merge(
        climate,
        (
            no3_deposition = fill(0.01f0, 3, 1),
            nh4_deposition = fill(0.02f0, 3, 1),
        ),
    )
    run_simulation!(simulation, deposition_climate; spinup = false)

    @test simulation.simulated_days == 3
    @test all(isfinite, simulation.daily_weather.no3_deposition)
    @test all(isfinite, simulation.daily_weather.nh4_deposition)

    diagnosed = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 3, fertilizer = :no,
        nitrogen_limit_vcmax = true, diagnostics = true,
    )
    run_simulation!(diagnosed, deposition_climate; spinup = false)
    balance = diagnosed.nitrogen_balance
    @test balance.atmospheric_deposition[:, 1] == fill(0.03f0, 3)
    @test maximum(abs, balance.residual[:, 1]) < 1.0f-4
end

@testset "Daily global CO₂ remains aligned across climate blocks" begin
    initial, _ = simulation_api_fixture(Float32)
    simulation = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 4, fertilizer = :no,
    )
    first_block = merge(
        climate_block(Float32, 2, 15, 1),
        (co2 = Float32[400, 401], co2_daily = true),
    )
    second_block = merge(
        climate_block(Float32, 2, 15, 1),
        (co2 = Float32[402, 403], co2_daily = true),
    )
    run_simulation!(simulation, first_block; spinup = false)
    run_simulation!(simulation, second_block; spinup = false)
    @test simulation.simulated_days == 4
    @test simulation.daily_weather.annual_co2[1] == 40.3f0
end

@testset "High-level C4 simulation uses the shared daily driver" begin
    initial, climate = simulation_api_fixture(Float32)
    simulation = initialize_simulation(
        cft4, initial;
        indices = [1], T = Float32, days = 3, fertilizer = :no,
    )
    run_simulation!(simulation, climate; spinup = false)
    @test simulation.simulated_days == 3
    @test size(simulation.output.crop.npp) == (3, 1)
    @test all(isfinite, simulation.output.crop.npp)
end

@testset "Daily ecosystem flux outputs close to their process fluxes" begin
    initial, climate = simulation_api_fixture(Float32)
    simulation = initialize_simulation(
        cft1, initial; indices = [1], T = Float32, days = 3, fertilizer = :no,
    )
    run_simulation!(simulation, climate; spinup = false)
    crop = simulation.state
    @test size(simulation.output.soil.ecosystem_respiration) == (3, 1)
    @test size(simulation.output.soil.heterotrophic_respiration) == (3, 1)
    @test size(simulation.output.soil.evapotranspiration) == (3, 1)
    @test simulation.output.soil.ecosystem_respiration[end, 1] ≈
          NeuralCrop.crop_fluxes(crop).carbon.respiration[1] +
          NeuralCrop.crop_fluxes(crop).carbon.leaf_respiration[1] +
          NeuralCrop.soil_carbon_fluxes(crop).heterotrophic_respiration[1]
    expected_et = NeuralCrop.crop_fluxes(crop).water.interception[1] +
        NeuralCrop.soil_surface_litter_fluxes(crop).evaporation[1] +
        sum(NeuralCrop.crop_fluxes(crop).water.transpiration_layer[:, 1]) +
        sum(NeuralCrop.soil_water_fluxes(crop).evaporation[:, 1])
    @test simulation.output.soil.evapotranspiration[end, 1] ≈ expected_et
    reco_spec, reco_frequency = output_variable_spec(:soil, :ecosystem_respiration)
    @test reco_spec.units == "gC m-2 day-1"
    @test reco_frequency === :daily

    chunks = OutputChunk[]
    streamed = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 3, diagnostics = false, fertilizer = :no,
    )
    stream = OutputStream(
        [OutputVariable(:soil, :evapotranspiration)];
        writer = chunk -> push!(chunks, chunk), cell_ids = Int32[1],
    )
    run_simulation!(streamed, [climate]; spinup = false, output_stream = stream)
    @test only(chunks).values[:soil_evapotranspiration] ==
          simulation.output.soil.evapotranspiration
end

@testset "Annual management is applied at sowing" begin
    initial, climate = simulation_api_fixture(Float32)
    simulation = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 1, fertilizer = :no,
    )
    management = (
        sdate = Int32[1],
        phu = Float32[-700],
        manure = Float32[2],
        fertilizer = Float32[3],
        residuefrac = Float32[0.5],
    )
    run_simulation!(simulation, climate; end_day = 1, spinup = false, management)
    @test simulation.state.inputs.crop.phenology.phu == Float32[700]
    @test simulation.state.inputs.crop.phenology.winter_type == Bool[true]
    @test simulation.managed_land.manure == Float32[2]
    @test simulation.managed_land.fertilizer == Float32[3]
    @test simulation.managed_land.residue_fraction == Float32[0.5]
end

@testset "One-day transition matches the range runner" begin
    initial, climate = simulation_api_fixture(Float32)
    baseline = initialize_simulation(
        cft1, initial; indices = [1], T = Float32, days = 3, fertilizer = :no,
    )
    stepped = initialize_simulation(
        cft1, initial; indices = [1], T = Float32, days = 3, fertilizer = :no,
    )
    run_simulation!(baseline, climate; spinup = false)
    for day in 1:3
        transition_day!(stepped, climate; climate_day = day)
    end
    @test stepped.simulated_days == 3
    @test stepped.output.crop.npp == baseline.output.crop.npp
    @test stepped.output.crop.biomass == baseline.output.crop.biomass
    @test stepped.state.prognostic.soil.water.storage ==
        baseline.state.prognostic.soil.water.storage
end

@testset "Runtime contracts describe active state and outputs" begin
    initial, _ = simulation_api_fixture(Float32)
    simulation = initialize_simulation(
        cft1, initial; indices = [1], T = Float32, days = 1, fertilizer = :no,
    )
    @test simulation.config isa SimulationConfiguration
    @test simulation.config.fertilizer === :no
    @test simulation.config.irrigation === false
    @test architecture_name(simulation.config.execution) == :cpu
    @test float_type(simulation.config.execution) === Float32
    @test simulation.config.execution.domain.cell_ids == Int32[1]
    @test !isempty(state_schema(simulation.state))
    @test validate_state_schema(simulation.state, 1) === simulation.state
    spec, frequency = output_variable_spec(:crop, :npp)
    @test spec.units == "gC m-2 day-1"
    @test frequency == :daily
end

@testset "Backend-neutral climate blocks use the simulation precision" begin
    initial, climate = simulation_api_fixture(Float32)
    simulation = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float64, days = 3, fertilizer = :no,
    )
    forcing = (
        temp = climate.temp,
        prec = climate.prec,
        sw = climate.swdown,
        lw = climate.lwnet,
        co2 = fill(400.0f0, 3),
        co2_daily = true,
        backend_neutral = true,
    )
    run_simulation!(simulation, [forcing]; spinup = false)
    @test simulation.simulated_days == 3
    @test eltype(simulation.daily_weather.temp) == Float64
    @test simulation.daily_weather.temp == [15.0]
end

@testset "Backend-neutral initial data need no source indices" begin
    initial, climate = simulation_api_fixture(Float32)
    backend_neutral = merge(initial, (backend_neutral = true,))
    simulation = initialize_simulation(
        cft1, backend_neutral;
        T = Float64, days = 3, fertilizer = :no,
    )
    forcing = (
        temp = climate.temp,
        prec = climate.prec,
        sw = climate.swdown,
        lw = climate.lwnet,
        co2 = fill(400.0f0, 3),
        co2_daily = true,
        backend_neutral = true,
    )
    run_simulation!(simulation, [forcing]; spinup = false)
    @test simulation.config.indices === nothing
    @test simulation.simulated_days == 3
    @test eltype(simulation.managed_land.latitude) == Float64
end

@testset "Backend-neutral initial data honor an optional selection" begin
    initial, _ = simulation_api_fixture(Float32)
    backend_neutral = merge(initial, (backend_neutral = true,))
    simulation = initialize_simulation(
        cft1, backend_neutral;
        indices = [1], T = Float32, days = 3, fertilizer = :no,
    )
    @test length(simulation.managed_land.latitude) == 1
    @test simulation.config.indices == [1]
end

@testset "Annual CO₂ forcing length is validated before kernel launch" begin
    initial, _ = simulation_api_fixture(Float32)
    simulation = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 366, fertilizer = :yes,
    )
    incomplete = merge(climate_block(Float32, 366, 15, 1), (co2 = Float32[400],))
    @test_throws DimensionMismatch run_simulation!(
        simulation, incomplete; spinup = false,
    )
end

@testset "High-level crop simulation API" begin
    initial, climate = simulation_api_fixture(Float32)
    simulation = initialize_simulation(
        cft1, initial;
        indices = [1],
        T = Float64,
        device = identity,
        days = 3,
        fertilizer = :yes,
    )

    @test simulation.output === simulation.state.output
    @test simulation.config.with_tillage
    @test simulation.processes.crop === simulation.cft
    @test simulation.processes.global_parameters === simulation.model_parameters
    crop_lifecycle_ids = runtime_array_ids((
        simulation.state.prognostic.crop,
        simulation.state.fluxes.crop,
        simulation.state.auxiliary.crop,
        simulation.state.inputs.crop,
        simulation.state.events.crop,
        simulation.state.workspace.crop,
    ))
    soil_lifecycle_ids = runtime_array_ids((
        simulation.state.prognostic.soil,
        simulation.state.fluxes.soil,
        simulation.state.auxiliary.soil,
        simulation.state.inputs.soil,
        simulation.state.workspace.soil,
    ))
    @test length(crop_lifecycle_ids) == length(unique(crop_lifecycle_ids))
    @test !isempty(crop_lifecycle_ids)
    @test length(soil_lifecycle_ids) == length(unique(soil_lifecycle_ids))
    @test !isempty(soil_lifecycle_ids)
    @test eltype(simulation.state.prognostic.crop.canopy.lai) == Float64
    @test eltype(simulation.water_balance.residual) == Float64

    returned = run_simulation!(simulation, climate; spinup = false)
    @test returned === simulation
    @test simulation.simulated_days == 3
    @test size(simulation.output.crop.npp) == (3, 1)
    @test all(isfinite, simulation.output.crop.npp)

    summary = simulation_summary(simulation)
    @test summary.precision == Float64
    @test summary.cells == 1
    @test summary.simulated_days == 3
    @test isfinite(summary.crop.cumulative_npp)
    @test_throws ArgumentError run_simulation!(simulation, climate; spinup = false)
end

@testset "Multiple climate blocks preserve the continuous daily timeline" begin
    initial, _ = simulation_api_fixture(Float32)
    first_block = climate_block(Float32, 2, 15, 1)
    second_block = climate_block(Float32, 2, 17, 3)
    continuous = (
        temp_spinup = first_block.temp_spinup,
        temp = vcat(first_block.temp, second_block.temp),
        prec = vcat(first_block.prec, second_block.prec),
        swdown = vcat(first_block.swdown, second_block.swdown),
        lwnet = vcat(first_block.lwnet, second_block.lwnet),
        windspeed = vcat(first_block.windspeed, second_block.windspeed),
        co2 = Float32[400],
    )

    create() = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float64, days = 4, fertilizer = :yes,
    )
    chunked = create()
    single = create()
    run_simulation!(chunked, [first_block, second_block]; spinup = false)
    run_simulation!(single, continuous; spinup = false)

    @test chunked.simulated_days == 4
    @test chunked.output.crop.npp ≈ single.output.crop.npp
    @test chunked.output.calendar.sowing_event == single.output.calendar.sowing_event
    @test findall(!iszero, vec(chunked.output.calendar.sowing_event)) == [1]
    @test chunked.state.prognostic.crop.carbon.leaf ≈ single.state.prognostic.crop.carbon.leaf
    @test chunked.state.prognostic.crop.carbon.root ≈ single.state.prognostic.crop.carbon.root
    @test chunked.state.prognostic.crop.carbon.pool ≈ single.state.prognostic.crop.carbon.pool
    @test chunked.state.prognostic.crop.carbon.storage ≈ single.state.prognostic.crop.carbon.storage
    @test chunked.state.prognostic.soil.water.storage ≈ single.state.prognostic.soil.water.storage
    @test chunked.water_balance.precipitation == reshape(Float64[1, 1, 3, 3], 4, 1)

    mktempdir() do directory
        first_path = joinpath(directory, "climate_1.jld2")
        second_path = joinpath(directory, "climate_2.jld2")
        JLD2.jldsave(first_path; climate = first_block)
        JLD2.jldsave(second_path; climate = second_block)
        from_files = create()
        run_simulation!(from_files, [first_path, second_path]; spinup = false)
        @test from_files.output.crop.npp ≈ single.output.crop.npp
        @test from_files.state.prognostic.soil.water.storage ≈ single.state.prognostic.soil.water.storage
    end
end

@testset "Checkpoint restore preserves continuous simulation" begin
    initial, _ = simulation_api_fixture(Float32)
    first_block = climate_block(Float32, 2, 15, 1)
    second_block = climate_block(Float32, 2, 17, 3)
    continuous = (
        temp_spinup = first_block.temp_spinup,
        temp = vcat(first_block.temp, second_block.temp),
        prec = vcat(first_block.prec, second_block.prec),
        swdown = vcat(first_block.swdown, second_block.swdown),
        lwnet = vcat(first_block.lwnet, second_block.lwnet),
        windspeed = vcat(first_block.windspeed, second_block.windspeed),
        co2 = Float32[400],
    )
    create(days = 4; cell_ids = [101]) = initialize_simulation(
        cft1, initial;
        indices = [1], cell_ids, T = Float32, days = days, fertilizer = :yes,
    )

    reference = create()
    run_simulation!(reference, continuous; spinup = false)

    interrupted = create()
    run_simulation!(interrupted, first_block; spinup = false)
    mktempdir() do directory
        path = joinpath(directory, "crop_checkpoint.jld2")
        @test save_checkpoint(path, interrupted) == path

        restored = create()
        @test restore_checkpoint!(restored, path) === restored
        @test restored.simulated_days == 2
        @test restored.output.crop.npp == interrupted.output.crop.npp
        @test restored.state.prognostic.crop.carbon.leaf == interrupted.state.prognostic.crop.carbon.leaf
        @test restored.state.prognostic.soil.water.storage == interrupted.state.prognostic.soil.water.storage
        @test restored.climbuf.atemp == interrupted.climbuf.atemp

        run_simulation!(restored, second_block; spinup = false)
        @test restored.simulated_days == reference.simulated_days
        @test restored.output.crop.npp == reference.output.crop.npp
        @test restored.output.crop.water_deficit ==
            reference.output.crop.water_deficit
        @test restored.state.prognostic.crop.carbon.leaf == reference.state.prognostic.crop.carbon.leaf
        @test restored.state.prognostic.crop.nitrogen.total == reference.state.prognostic.crop.nitrogen.total
        @test restored.state.prognostic.soil.water.storage == reference.state.prognostic.soil.water.storage
        @test restored.state.prognostic.soil.carbon.fast == reference.state.prognostic.soil.carbon.fast
        @test restored.state.prognostic.soil.nitrogen.nitrate == reference.state.prognostic.soil.nitrogen.nitrate
        @test restored.water_balance.residual == reference.water_balance.residual
        @test restored.nitrogen_balance.residual == reference.nitrogen_balance.residual
        @test restored.carbon_balance.residual == reference.carbon_balance.residual
        @test restored.thermal_balance.energy_residual ==
            reference.thermal_balance.energy_residual

        incompatible = create(5)
        @test_throws ArgumentError restore_checkpoint!(incompatible, path)

        wrong_domain = create(; cell_ids = [202])
        @test_throws ArgumentError restore_checkpoint!(wrong_domain, path)

        checkpoint = JLD2.load(path, "checkpoint")
        incompatible_metadata = merge(
            checkpoint.metadata,
            (parameter_fingerprint = "incompatible-parameters",),
        )
        JLD2.jldsave(path; checkpoint = merge(
            checkpoint,
            (metadata = incompatible_metadata,),
        ))
        @test_throws ArgumentError restore_checkpoint!(create(), path)
    end
end
