using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

const C3_E2E_DAYS = 730
const C3_E2E_CELLS = 2

to_backend(::typeof(identity), values) = copy(values)
to_backend(device, values) = device(values)

function c3_e2e_initial_data(device)
    cells = C3_E2E_CELLS
    layers = 5

    crop = (
        sdate = to_backend(device, Int32[100, 105]),
        phu = to_backend(device, Float32[543, 620]),
        manure = to_backend(device, zeros(Float32, cells)),
        fertilizer = to_backend(device, Float32[24.55, 30]),
        residuefrac = to_backend(device, Float32[0.67, 0.55]),
    )
    soil_parameters = (
        ph = to_backend(device, Float32[6.5, 6.2]),
        w_sat = to_backend(device, repeat(Float32[0.404, 0.43]', layers, 1)),
        sand = to_backend(device, reshape(Float32[0.58, 0.46], 1, cells)),
        clay = to_backend(device, reshape(Float32[0.27, 0.22], 1, cells)),
        tdiff_0 = to_backend(device, Float32[0.78, 0.70]),
        tdiff_15 = to_backend(device, Float32[0.808, 0.75]),
        # Keep this host-side: init_soil transfers the immutable layer geometry.
        soildepth = Float32[200, 300, 500, 1000, 1000],
    )

    swc = Float32[57.41 62.0; 55.32 68.0; 126.13 140.0; 274.59 260.0; 285.71 275.0]
    litc = Float32[0.13 0.2; 187.5 165.0; 225.36 205.0]
    fastc = Float32[548.97 510.0; 368.27 350.0; 313.79 300.0; 377.55 360.0; 344.65 330.0]
    slowc = Float32[1218.62 1150.0; 753.33 720.0; 660.10 640.0; 792.63 760.0; 736.38 710.0]
    litn = Float32[0.0047 0.007; 6.47 5.8; 9.47 8.7]
    fastn = Float32[36.60 34.0; 24.55 23.0; 20.92 20.0; 25.17 24.0; 22.98 22.0]
    slown = Float32[81.24 77.0; 50.22 48.0; 44.01 42.5; 52.84 50.5; 49.09 47.0]
    u0 = (
        swc = to_backend(device, swc),
        litc = to_backend(device, litc),
        fastc = to_backend(device, fastc),
        slowc = to_backend(device, slowc),
        litn = to_backend(device, litn),
        fastn = to_backend(device, fastn),
        slown = to_backend(device, slown),
    )

    return (
        latitude = to_backend(device, Float32[45, 47]),
        soilparams = soil_parameters,
        ModelState = (crop = crop, u0 = u0),
    )
end

function c3_e2e_climate_host()
    days = C3_E2E_DAYS
    cells = C3_E2E_CELLS
    temperature = zeros(Float32, days, cells)
    precipitation = zeros(Float32, days, cells)
    shortwave = zeros(Float32, days, cells)
    longwave = zeros(Float32, days, cells)
    wind = zeros(Float32, days, cells)

    for cell in 1:cells, day in 1:days
        day_of_year = mod1(day, 365)
        phase = 2.0f0 * Float32(pi) * Float32(day_of_year - 109 - 5 * (cell - 1)) / 365.0f0
        seasonal = sin(phase)
        temperature[day, cell] = 9.0f0 + 11.0f0 * seasonal + 0.4f0 * (cell - 1)
        shortwave[day, cell] = max(20.0f0, 175.0f0 + 125.0f0 * seasonal)
        longwave[day, cell] = -45.0f0 + 12.0f0 * seasonal
        wind[day, cell] = 1.8f0 + 0.2f0 * (cell - 1)
        precipitation[day, cell] = if mod(day + 2 * cell, 11) == 0
            11.0f0
        elseif mod(day + cell, 5) == 0
            2.5f0
        else
            0.0f0
        end
    end

    return (
        temp = temperature,
        prec = precipitation,
        sw = shortwave,
        lw = longwave,
        wind = wind,
        co2 = Float32[370, 373],
    )
end

function c3_e2e_climate(device, host)
    return (
        temp = to_backend(device, host.temp),
        prec = to_backend(device, host.prec),
        sw = to_backend(device, host.sw),
        lw = to_backend(device, host.lw),
        wind = to_backend(device, host.wind),
        co2 = to_backend(device, host.co2),
    )
end

function initialize_c3_e2e_climbuf!(climbuf, host_climate, device)
    # Use the same deterministic one-year history on both backends without
    # adding another 365 kernel launches to this end-to-end regression.
    climbuf.atemp .= to_backend(device, host_climate.temp[1:365, :])
    climbuf.temp .= to_backend(device, host_climate.temp[335:365, :])
    climbuf.atemp_mean .= to_backend(
        device,
        vec(sum(host_climate.temp[1:365, :]; dims = 1) ./ 365.0f0),
    )
    climbuf.V_req .= 0.0f0
    return nothing
end

function run_c3_e2e(device, host_climate)
    initial_data = c3_e2e_initial_data(device)
    climbuf, crop, pet, soil, managed_land, weather, output = init_states!(
        cft1, initial_data, C3_E2E_CELLS, device,
    )
    climate = c3_e2e_climate(device, host_climate)
    initialize_c3_e2e_climbuf!(climbuf, host_climate, device)

    water = init_water_balance(C3_E2E_DAYS, C3_E2E_CELLS, device)
    nitrogen = init_nitrogen_balance(C3_E2E_DAYS, C3_E2E_CELLS, device)
    carbon = init_carbon_balance(C3_E2E_DAYS, C3_E2E_CELLS, device)
    thermal = init_thermal_balance(C3_E2E_DAYS, C3_E2E_CELLS, device)

    state = model_state(climbuf, crop, pet, soil, managed_land, weather, output)
    processes = ProcessModules(cft1, ModelParameters(Float32))
    daily_crop_C3!(
        1, C3_E2E_DAYS, processes, climate, state;
        irrigation = false,
        manure = false,
        fertilizer = :yes,
        nitrogen_limit_vcmax = false,
        water_balance = water,
        nitrogen_balance = nitrogen,
        carbon_balance = carbon,
        thermal_balance = thermal,
    )

    return (; state, climbuf, crop, pet, soil, managed_land, weather, output,
            water, nitrogen, carbon, thermal)
end

host_array(values) = Array(values)

function test_float_equivalence(gpu_values, cpu_values;
                                rtol = 1.0f-3, atol = 5.0f-4,
                                label = "unnamed")
    gpu_host = host_array(gpu_values)
    cpu_host = host_array(cpu_values)
    @test size(gpu_host) == size(cpu_host)
    @test all(isfinite, gpu_host)
    @test all(isfinite, cpu_host)
    # Array `isapprox` uses a norm, so thousands of individually acceptable
    # daily Float32 differences can fail after accumulation. Model equivalence
    # is a pointwise requirement: every day/cell must satisfy the tolerance.
    matches = isapprox.(gpu_host, cpu_host; rtol = rtol, atol = atol)
    if !all(matches)
        first_bad = findfirst(.!matches)
        gpu_value = gpu_host[first_bad]
        cpu_value = cpu_host[first_bad]
        absolute_error = abs(gpu_value - cpu_value)
        relative_error = absolute_error / max(abs(cpu_value), atol)
        @info "CPU/GPU first pointwise mismatch" label first_bad gpu_value cpu_value absolute_error relative_error rtol atol
    end
    @test all(matches)
    return nothing
end

function test_balance_equivalence(gpu_balance, cpu_balance;
                                  group, rtol, atol,
                                  field_rtol = NamedTuple(),
                                  field_atol = NamedTuple(),
                                  skip_fields = ())
    for field in fieldnames(typeof(cpu_balance))
        field in skip_fields && continue
        comparison_rtol = hasproperty(field_rtol, field) ?
            getproperty(field_rtol, field) : rtol
        comparison_atol = hasproperty(field_atol, field) ?
            getproperty(field_atol, field) : atol
        test_float_equivalence(
            getproperty(gpu_balance, field),
            getproperty(cpu_balance, field);
            rtol = comparison_rtol,
            atol = comparison_atol,
            label = "$group.$field",
        )
    end
    return nothing
end

function test_thermal_closure(balance; label)
    column_energy = host_array(balance.column_energy)
    surface_energy = host_array(balance.surface_energy_flux)
    energy_residual = host_array(balance.energy_residual)
    energy_scale = max.(abs.(column_energy), abs.(surface_energy), 1.0f0)
    # The five-layer column is O(1e8) J m-2 in Float32. Its daily closure is
    # therefore tested against a small multiple of the ledger's representable
    # spacing, rather than comparing two cancellation-sensitive residuals.
    energy_tolerance = 16.0f0 .* eps(Float32) .* energy_scale .+ 1.0f0
    maximum_absolute_residual = maximum(abs, energy_residual)
    maximum_relative_residual = maximum(abs.(energy_residual) ./ energy_scale)
    @info "Thermal energy closure" label maximum_absolute_residual maximum_relative_residual
    @test all(isfinite, energy_residual)
    @test all(abs.(energy_residual) .<= energy_tolerance)

    percolation_residual = host_array(balance.percolation_energy_residual)
    boundary_scale = max.(
        abs.(host_array(balance.rain_energy_input)) .+
        abs.(host_array(balance.snowmelt_energy_input)) .+
        abs.(host_array(balance.lateral_runoff_energy_output)) .+
        abs.(host_array(balance.bottom_drainage_energy_output)),
        1.0f0,
    )
    percolation_tolerance = 5.0f-6 .* boundary_scale .+ 2.0f0
    maximum_absolute_residual = maximum(abs, percolation_residual)
    maximum_relative_residual = maximum(abs.(percolation_residual) ./ boundary_scale)
    @info "Percolation energy closure" label maximum_absolute_residual maximum_relative_residual
    @test all(isfinite, percolation_residual)
    @test all(abs.(percolation_residual) .<= percolation_tolerance)
    return nothing
end

function event_days(values, cell)
    host = host_array(values)
    return findall(!iszero, view(host, :, cell))
end

function log_first_crop_divergence(cpu, gpu;
                                   rtol = 1.0f-5, atol = 1.0f-6)
    for cell in 1:C3_E2E_CELLS
        cpu_sowing = event_days(cpu.output.calendar.sowing_event, cell)
        gpu_sowing = event_days(gpu.output.calendar.sowing_event, cell)
        cpu_harvest = event_days(cpu.output.calendar.harvest_event, cell)
        gpu_harvest = event_days(gpu.output.calendar.harvest_event, cell)
        @info "CPU/GPU crop event days" cell cpu_sowing gpu_sowing cpu_harvest gpu_harvest
    end

    # Follow the daily causal chain from phenology/photosynthesis to crop state.
    for field in (:fphu, :lambda, :gpp, :respiration, :npp, :biomass, :lai)
        cpu_values = host_array(getproperty(cpu.output.crop, field))
        gpu_values = host_array(getproperty(gpu.output.crop, field))
        matches = isapprox.(gpu_values, cpu_values; rtol = rtol, atol = atol)
        first_bad = findfirst(.!matches)
        first_bad === nothing && continue
        gpu_value = gpu_values[first_bad]
        cpu_value = cpu_values[first_bad]
        absolute_error = abs(gpu_value - cpu_value)
        relative_error = absolute_error / max(abs(cpu_value), atol)
        @info "CPU/GPU key crop divergence" field first_bad gpu_value cpu_value absolute_error relative_error rtol atol
    end
    return nothing
end

function test_exact_equivalence(gpu_values, cpu_values; label)
    gpu_host = host_array(gpu_values)
    cpu_host = host_array(cpu_values)
    matches = gpu_host .== cpu_host
    if !all(matches)
        first_bad = findfirst(.!matches)
        gpu_value = gpu_host[first_bad]
        cpu_value = cpu_host[first_bad]
        @info "CPU/GPU first exact mismatch" label first_bad gpu_value cpu_value
    end
    @test all(matches)
    return nothing
end

"""Recursively compare every numerical leaf in the canonical lifecycle tree."""
function test_lifecycle_state_equivalence(gpu_value, cpu_value; label)
    # These are cancellation-sensitive diagnostic residuals, not prognostic
    # state. Each backend is instead checked against its own energy ledger by
    # `test_thermal_closure` below.
    if label in (
        "state.fluxes.soil.thermal.energy_residual",
        "state.fluxes.soil.thermal.percolation_energy_residual",
    )
        return nothing
    end

    if gpu_value isa AbstractArray
        if eltype(gpu_value) <: AbstractFloat
            test_float_equivalence(gpu_value, cpu_value; label)
        else
            test_exact_equivalence(gpu_value, cpu_value; label)
        end
        return nothing
    end

    if gpu_value isa NamedTuple || isstructtype(typeof(gpu_value))
        @test propertynames(gpu_value) == propertynames(cpu_value)
        for field in propertynames(cpu_value)
            test_lifecycle_state_equivalence(
                getproperty(gpu_value, field), getproperty(cpu_value, field);
                label = "$label.$field",
            )
        end
        return nothing
    end

    @test gpu_value == cpu_value
    return nothing
end

@testset "CUDA C3 rainfed wheat 365/730-day end-to-end equivalence" begin
    host_climate = c3_e2e_climate_host()
    cpu = run_c3_e2e(identity, host_climate)
    gpu = run_c3_e2e(CuArray, host_climate)
    synchronize()

    @test size(cpu.output.crop.npp) == (730, C3_E2E_CELLS)
    @test size(gpu.output.crop.npp) == (730, C3_E2E_CELLS)
    @test maximum(cpu.output.crop.npp) > 0.0f0
    @test sum(cpu.output.calendar.sowing_event) > 0
    @test sum(cpu.output.calendar.harvest_event) > 0

    test_lifecycle_state_equivalence(gpu.state, cpu.state; label = "state")
    log_first_crop_divergence(cpu, gpu)

    daily_float_fields = (
        :gpp, :npp, :lambda, :potential_vcmax, :vcmax,
        :nitrogen_limitation, :respiration, :biomass, :lai,
        :storage_carbon, :fphu, :water_deficit,
    )
    for field in daily_float_fields
        cpu_values = getproperty(cpu.output.crop, field)
        gpu_values = getproperty(gpu.output.crop, field)
        # Explicit one-year checkpoint and complete two-year trajectory.
        test_float_equivalence(
            gpu_values[1:365, :], cpu_values[1:365, :];
            label = "output.crop.$field[1:365]",
        )
        test_float_equivalence(
            gpu_values, cpu_values; label = "output.crop.$field[1:730]",
        )
    end
    test_float_equivalence(
        gpu.output.crop.yield, cpu.output.crop.yield;
        label = "output.crop.yield",
    )

    for field in (:growing_mask,)
        test_exact_equivalence(
            getproperty(gpu.output.crop, field),
            getproperty(cpu.output.crop, field);
            label = "output.crop.$field",
        )
    end
    for field in (
        :harvesting_mask, :harvesting_year, :harvest_date,
        :sowing_event, :harvest_event,
    )
        test_exact_equivalence(
            getproperty(gpu.output.calendar, field),
            getproperty(cpu.output.calendar, field);
            label = "output.calendar.$field",
        )
    end

    crop_state_fields = (
        ("crop.state.carbon.leaf", cpu.crop.state.carbon.leaf, gpu.crop.state.carbon.leaf),
        ("crop.state.carbon.root", cpu.crop.state.carbon.root, gpu.crop.state.carbon.root),
        ("crop.state.carbon.pool", cpu.crop.state.carbon.pool, gpu.crop.state.carbon.pool),
        ("crop.state.carbon.storage", cpu.crop.state.carbon.storage, gpu.crop.state.carbon.storage),
        ("crop.state.nitrogen.total", cpu.crop.state.nitrogen.total, gpu.crop.state.nitrogen.total),
        ("crop.state.nitrogen.leaf", cpu.crop.state.nitrogen.leaf, gpu.crop.state.nitrogen.leaf),
        ("crop.state.nitrogen.root", cpu.crop.state.nitrogen.root, gpu.crop.state.nitrogen.root),
        ("crop.state.nitrogen.pool", cpu.crop.state.nitrogen.pool, gpu.crop.state.nitrogen.pool),
        ("crop.state.nitrogen.storage", cpu.crop.state.nitrogen.storage, gpu.crop.state.nitrogen.storage),
        ("crop.state.canopy.lai", cpu.crop.state.canopy.lai, gpu.crop.state.canopy.lai),
        ("crop.auxiliary.phenology.fphu", cpu.crop.auxiliary.phenology.fphu, gpu.crop.auxiliary.phenology.fphu),
    )
    for (label, cpu_values, gpu_values) in crop_state_fields
        test_float_equivalence(gpu_values, cpu_values; label = label)
    end

    soil_state_fields = (
        ("soil.water.storage", cpu.soil.water.storage, gpu.soil.water.storage),
        ("soil.water.ice_storage", cpu.soil.water.ice_storage, gpu.soil.water.ice_storage),
        ("soil.thermal.temperature", cpu.soil.thermal.temperature, gpu.soil.thermal.temperature),
        ("soil.thermal.enthalpy", cpu.soil.thermal.enthalpy, gpu.soil.thermal.enthalpy),
        ("soil.carbon.litter", cpu.soil.carbon.litter, gpu.soil.carbon.litter),
        ("soil.carbon.fast", cpu.soil.carbon.fast, gpu.soil.carbon.fast),
        ("soil.carbon.slow", cpu.soil.carbon.slow, gpu.soil.carbon.slow),
        ("soil.nitrogen.litter", cpu.soil.nitrogen.litter, gpu.soil.nitrogen.litter),
        ("soil.nitrogen.fast", cpu.soil.nitrogen.fast, gpu.soil.nitrogen.fast),
        ("soil.nitrogen.slow", cpu.soil.nitrogen.slow, gpu.soil.nitrogen.slow),
        ("soil.nitrogen.nitrate", cpu.soil.nitrogen.nitrate, gpu.soil.nitrogen.nitrate),
        ("soil.nitrogen.ammonium", cpu.soil.nitrogen.ammonium, gpu.soil.nitrogen.ammonium),
        ("soil.management.tillage_density_factor",
         cpu.soil.management.tillage_density_factor,
         gpu.soil.management.tillage_density_factor),
    )
    for (label, cpu_values, gpu_values) in soil_state_fields
        test_float_equivalence(gpu_values, cpu_values; label = label)
    end
    @test host_array(gpu.soil.thermal.initialized) ==
        host_array(cpu.soil.thermal.initialized)

    test_balance_equivalence(
        gpu.water, cpu.water; group = "water", rtol = 1.0f-3, atol = 5.0f-4,
    )
    test_balance_equivalence(
        gpu.nitrogen, cpu.nitrogen;
        group = "nitrogen", rtol = 2.0f-3, atol = 1.0f-3,
    )
    test_balance_equivalence(
        gpu.carbon, cpu.carbon;
        group = "carbon", rtol = 2.0f-3, atol = 1.0f-3,
        # The residual subtracts stocks of thousands of g C m-2. A 2 mg C m-2
        # CPU/GPU tolerance covers only final-bit Float32 cancellation.
        field_atol = (residual = 2.0f-3,),
    )
    test_balance_equivalence(
        gpu.thermal, cpu.thermal;
        group = "thermal", rtol = 2.0f-3, atol = 1.0f-1,
        # A single drainage pulse can cross a threshold on adjacent days as
        # CPU/GPU Float32 trajectories round differently; compare its
        # accumulated energy below rather than its daily pulse.
        field_atol = (untracked_water_energy_flux = 8.0f0,),
        # Residuals are conservation tests, not trajectory state: validate each
        # backend against its own physical ledger below.
        skip_fields = (
            :energy_residual,
            :percolation_energy_residual,
            :bottom_drainage_energy_output,
        ),
    )
    test_float_equivalence(
        sum(host_array(gpu.thermal.bottom_drainage_energy_output); dims = 1),
        sum(host_array(cpu.thermal.bottom_drainage_energy_output); dims = 1);
        label = "thermal.accumulated_bottom_drainage_energy_output",
        rtol = 5.0f-3,
        atol = 1.0f-1,
    )
    test_thermal_closure(cpu.thermal; label = "CPU")
    test_thermal_closure(gpu.thermal; label = "GPU")
end
