"""
    CropSimulation

High-level simulation container. It groups model state, diagnostics, numerical
precision, backend, and run options. Numerical variables are accessed through
the lifecycle groups in `simulation.state`.
"""
mutable struct CropSimulation{P, S, D, C}
    processes::P        # Process choices and immutable parameters; contains no runtime arrays.
    state::S            # Numerical arrays grouped uniformly by lifecycle.
    diagnostics::D      # Optional daily C/N/water/energy balance ledgers.
    config::C           # Backend, precision, domain-selection, and run-option metadata.
    simulated_days::Int # Number of completed daily time steps in this simulation.
end

function Base.show(io::IO, simulation::CropSimulation)
    diagnostics_enabled = any(!isnothing, values(simulation.diagnostics))
    print(
        io,
        "CropSimulation(", CFT_NAMES[simulation.cft.name], ", ",
        architecture_name(simulation.config.execution), ", ", simulation.config.T,
        ", cells=", length(simulation.config.execution.domain.indices),
        ", days=", simulation.simulated_days, "/", simulation.config.days,
        ", diagnostics=", diagnostics_enabled, ")",
    )
end

const _SIMULATION_STATE_PROPERTIES = (
    :climbuf, :pet, :managed_land, :daily_weather, :output,
)
const _SIMULATION_DIAGNOSTIC_PROPERTIES = (
    :water_balance, :nitrogen_balance, :carbon_balance, :thermal_balance,
)

function Base.getproperty(simulation::CropSimulation, name::Symbol)
    if name === :cft
        return getfield(simulation, :processes).crop
    elseif name === :model_parameters
        return getfield(simulation, :processes).global_parameters
    elseif name === :climbuf
        return getfield(simulation, :state).prognostic.climate
    elseif name === :pet
        return getfield(simulation, :state).auxiliary.pet
    elseif name === :managed_land
        return getfield(simulation, :state).inputs.management
    elseif name === :daily_weather
        return getfield(simulation, :state).inputs.weather
    elseif name === :output
        return getfield(simulation, :state).output
    elseif name in _SIMULATION_DIAGNOSTIC_PROPERTIES
        return getproperty(getfield(simulation, :diagnostics), name)
    end
    return getfield(simulation, name)
end

function Base.propertynames(::CropSimulation, private::Bool = false)
    public = (
        :processes, :cft, :model_parameters, :config, :simulated_days,
        :state, _SIMULATION_STATE_PROPERTIES..., _SIMULATION_DIAGNOSTIC_PROPERTIES...,
    )
    return private ? (public..., :diagnostics) : public
end

function _has_c_shift_restart(data::NamedTuple)
    return (hasproperty(data, :initial_state) &&
            hasproperty(data.initial_state, :c_shift_fast) &&
            hasproperty(data.initial_state, :c_shift_slow)) ||
        (hasproperty(data, :c_shift_fast) && hasproperty(data, :c_shift_slow))
end

function _prepare_initial_data(data::NamedTuple, indices, device, T)
    if hasproperty(data, :ModelState) && hasproperty(data, :soilparams)
        return data
    end
    if hasproperty(data, :backend_neutral) && data.backend_neutral
        restore_c_shift = _has_c_shift_restart(data)
        selected_indices = indices === nothing ? collect(1:length(data.latitude)) :
            collect(Int, indices)
        return InitialDataLoader(
            data, selected_indices, device;
            T = T,
            load_c_shift_restart = restore_c_shift,
        )
    end
    indices === nothing && throw(ArgumentError(
        "indices are required when initialize_simulation receives raw initial data",
    ))
    restore_c_shift = _has_c_shift_restart(data)
    return InitialDataLoader(
        data, collect(Int, indices), device;
        T = T,
        load_c_shift_restart = restore_c_shift,
    )
end

"""
    initialize_simulation(cft, initial_data; kwargs...)

Create a precision- and backend-consistent crop simulation. `initial_data` may
be either raw input data accepted by `InitialDataLoader` or its normalized
output. Set `diagnostics=false` to avoid allocating daily balance ledgers.
`fertilizer` follows LPJmL's `:no`, `:yes`, and `:auto` modes; `manure` is an
independent prescribed-input switch. By default,
`c_shift_initialization=:auto` restores supplied `c_shift_fast`/`c_shift_slow`
arrays and otherwise uses the model default vertical distribution.
`sowing_mode=:prescribed_sdate` uses management sowing dates. The optional
`:dynamic_sdate` mode derives a climate-calendar candidate with the LPJmL CFT
sowing method (winter wheat: temperature/winter type; maize:
temperature/precipitation; rice and soybean: precipitation), then applies its
bounded anomaly relative to the first valid climate calendar to the prescribed
spatial date.
"""
function initialize_simulation(
    cft::CFTParameters,
    initial_data::NamedTuple;
    indices = nothing,
    cell_ids = nothing,
    device = identity,
    T::Type{<:AbstractFloat} = Float32,
    days::Integer,
    diagnostics::Bool = true,
    irrigation::Bool = false,
    manure::Bool = false,
    fertilizer = :auto,
    with_tillage::Bool = true,
    crop_resp_fix::Bool = false,
    nitrogen_limit_vcmax::Bool = false,
    freeze_vernalization_requirement::Bool = false,
    sowing_mode::Symbol = :prescribed_sdate,
    mineral_nitrogen_initialization::Symbol = :from_slow_organic_nitrogen,
    c_shift_initialization::Symbol = :auto,
)
    days > 0 || throw(ArgumentError("days must be positive"))
    fertilizer = fertilizer_mode(fertilizer)
    sowing_mode in (:prescribed_sdate, :dynamic_sdate) || throw(ArgumentError(
        "sowing_mode must be :prescribed_sdate or :dynamic_sdate",
    ))
    sowing_mode !== :dynamic_sdate || cft.name in (Int32(1), Int32(2), Int32(3), Int32(9)) ||
        throw(ArgumentError(
            ":dynamic_sdate is currently implemented for wheat, rice, maize, and soybean only",
        ))
    prepared = _prepare_initial_data(initial_data, indices, device, T)
    cells = length(prepared.latitude)
    resolved_c_shift_initialization = if c_shift_initialization === :auto
        hasproperty(prepared.ModelState, :c_shift_fast) ? :restart : :default
    else
        c_shift_initialization
    end
    states = init_states!(
        cft, prepared, cells, device;
        T = T,
        mineral_nitrogen_initialization = mineral_nitrogen_initialization,
        c_shift_initialization = resolved_c_shift_initialization,
    )
    climbuf, crop, pet, soil, managed_land, daily_weather, output = states
    state = model_state(climbuf, crop, pet, soil, managed_land, daily_weather, output)
    active_indices = indices === nothing ? collect(1:cells) : collect(Int, indices)
    active_cell_ids = cell_ids === nothing ? active_indices : collect(Int, cell_ids)
    config = SimulationConfiguration(
        T, device, days, active_indices, active_cell_ids;
        indices,
        irrigation,
        manure,
        fertilizer,
        with_tillage,
        crop_resp_fix,
        nitrogen_limit_vcmax,
        freeze_vernalization_requirement,
        sowing_mode,
    )
    validate_state_schema(state, cells)
    balances = diagnostics ? (
        water_balance = init_water_balance(days, cells, device; T = T),
        nitrogen_balance = init_nitrogen_balance(days, cells, device; T = T),
        carbon_balance = init_carbon_balance(days, cells, device; T = T),
        thermal_balance = init_thermal_balance(days, cells, device; T = T),
    ) : (
        water_balance = nothing,
        nitrogen_balance = nothing,
        carbon_balance = nothing,
        thermal_balance = nothing,
    )
    processes = ProcessModules(convert_precision(T, cft), ModelParameters(T))
    return CropSimulation(
        processes, state, balances, config, 0,
    )
end

function _prepare_climate(simulation::CropSimulation, climate::NamedTuple)
    if hasproperty(climate, :backend_neutral) && climate.backend_neutral
        T = simulation.config.T
        prepared = (
            temp = T.(climate.temp),
            prec = T.(climate.prec),
            sw = T.(climate.sw),
            lw = T.(climate.lw),
            co2 = T.(climate.co2),
            co2_daily = true,
        )
        if hasproperty(climate, :wind)
            prepared = merge(prepared, (wind = T.(climate.wind),))
        end
        for name in (:no3_deposition, :nh4_deposition)
            hasproperty(climate, name) || continue
            values = T.(getproperty(climate, name))
            prepared = merge(
                prepared,
                (; name => _daily_deposition_forcing(values, size(prepared.temp, 1))),
            )
        end
        if hasproperty(climate, :temp_spinup)
            prepared = merge(prepared, (temp_spinup = T.(climate.temp_spinup),))
        end
        if hasproperty(climate, :prec_spinup)
            prepared = merge(prepared, (prec_spinup = T.(climate.prec_spinup),))
        end
        return _adapt_to_device(simulation.config.device, prepared)
    end
    if hasproperty(climate, :swdown) && hasproperty(climate, :lwnet)
        indices = simulation.config.indices
        indices === nothing && throw(ArgumentError(
            "raw climate data require indices in initialize_simulation",
        ))
        return ClimateDataLoader(
            climate, indices, simulation.config.device; T = simulation.config.T,
        )
    end
    return climate
end

function _prepare_annual_management!(simulation::CropSimulation, management::NamedTuple)
    required = (:sdate, :phu, :manure, :fertilizer, :residuefrac)
    all(hasproperty(management, name) for name in required) || throw(ArgumentError(
        "annual management requires fields $(join(string.(required), ", "))",
    ))
    cells = length(simulation.config.execution.domain.indices)
    all(length(getproperty(management, name)) == cells for name in required) ||
        throw(DimensionMismatch("every annual management field must contain $cells cells"))
    T = simulation.config.T
    device = simulation.config.device
    to_float(values) = device(T.(vec(values)))
    to_integer(values) = device(Int32.(round.(vec(values))))

    signed_phu = to_float(management.phu)
    prescribed_phu = abs.(signed_phu)
    prescribed_winter_type = signed_phu .< zero(T)
    all(isfinite, prescribed_phu) || throw(ArgumentError(
        "annual PHU contains non-finite values",
    ))
    calendar = crop_calendar_input(simulation.state)
    calendar.prescribed_sowing_date .= to_integer(management.sdate)
    if simulation.config.sowing_mode === :prescribed_sdate
        calendar.sowing_date .= calendar.prescribed_sowing_date
    end
    simulation.managed_land.manure .= to_float(management.manure)
    simulation.managed_land.fertilizer .= to_float(management.fertilizer)
    simulation.managed_land.residue_fraction .= to_float(management.residuefrac)
    return (; prescribed_phu, prescribed_winter_type)
end

function _transition_range!(
    simulation::CropSimulation,
    prepared_climate::NamedTuple,
    start_day::Integer,
    end_day::Integer,
    ; reuse_output::Bool = false,
    simulation_day_offset::Integer = simulation.simulated_days,
    selected_output::Union{Nothing, Set{Tuple{Symbol, Symbol}}} = nothing,
    prescribed_phu = nothing,
    prescribed_winter_type = nothing,
    neural_parameters = nothing,
    neural_layout = nothing,
)
    common = (
        irrigation = simulation.config.irrigation,
        manure = simulation.config.manure,
        fertilizer = simulation.config.fertilizer,
        with_tillage = simulation.config.with_tillage,
        crop_resp_fix = simulation.config.crop_resp_fix,
        nitrogen_limit_vcmax = simulation.config.nitrogen_limit_vcmax,
        update_vernalization_requirement = !simulation.config.freeze_vernalization_requirement,
        sowing_mode = simulation.config.sowing_mode,
        water_balance = simulation.water_balance,
        nitrogen_balance = simulation.nitrogen_balance,
        carbon_balance = simulation.carbon_balance,
        thermal_balance = simulation.thermal_balance,
        simulation_day_offset = simulation_day_offset,
        diagnostic_offset = simulation.simulated_days,
        reuse_output = reuse_output,
        selected_output = selected_output,
        prescribed_phu = prescribed_phu,
        prescribed_winter_type = prescribed_winter_type,
        neural_parameters = neural_parameters,
        neural_layout = neural_layout,
    )
    _daily_crop!(
        pathway_value(simulation.processes.pathway),
        start_day, end_day, simulation.processes, prepared_climate,
        simulation.state;
        common...,
    )
    return simulation
end

"""
    transition_day!(simulation, climate; climate_day=1)

Advance exactly one daily transition. This is the narrow numerical lifecycle
boundary used by the high-level runner and the future differentiable path.
"""
function transition_day!(
    simulation::CropSimulation,
    climate::NamedTuple;
    climate_day::Integer = 1,
    management::Union{Nothing, NamedTuple} = nothing,
    neural_parameters = nothing,
    neural_layout = nothing,
)
    simulation.simulated_days < simulation.config.days || throw(ArgumentError(
        "simulation already contains all $(simulation.config.days) configured days",
    ))
    prepared_climate = _prepare_climate(simulation, climate)
    1 <= climate_day <= size(prepared_climate.temp, 1) || throw(BoundsError(
        prepared_climate.temp, (climate_day, :),
    ))
    prescribed = isnothing(management) ?
        (prescribed_phu = nothing, prescribed_winter_type = nothing) :
        _prepare_annual_management!(simulation, management)
    _transition_range!(
        simulation, prepared_climate, climate_day, climate_day;
        simulation_day_offset = simulation.simulated_days + 1 - climate_day,
        neural_parameters,
        neural_layout,
        prescribed...,
    )
    simulation.simulated_days += 1
    return simulation
end

"""
    run_simulation!(simulation, climate; start_day=1, end_day=nothing, spinup=true)

Append one continuous climate block. File-local climate rows start at one for
every block, while phenology, outputs, and diagnostics continue from
`simulation.simulated_days`.
"""
function run_simulation!(
    simulation::CropSimulation,
    climate::NamedTuple;
    start_day::Integer = 1,
    end_day::Union{Nothing, Integer} = nothing,
    spinup::Bool = true,
    spinup_years::Integer = 1,
    reuse_output::Bool = false,
    selected_output::Union{Nothing, Set{Tuple{Symbol, Symbol}}} = nothing,
    management::Union{Nothing, NamedTuple} = nothing,
    neural_parameters = nothing,
    neural_layout = nothing,
)
    prepared_climate = _prepare_climate(simulation, climate)
    climate_days = size(prepared_climate.temp, 1)
    remaining_days = simulation.config.days - simulation.simulated_days
    remaining_days > 0 || throw(ArgumentError(
        "simulation already contains all $(simulation.config.days) configured days",
    ))
    local_end_day = end_day === nothing ?
        min(climate_days, start_day + remaining_days - 1) : Int(end_day)
    1 <= start_day <= local_end_day <= climate_days || throw(ArgumentError(
        "require 1 <= start_day <= end_day <= $climate_days climate rows",
    ))
    if ndims(prepared_climate.co2) == 1
        co2_daily = hasproperty(prepared_climate, :co2_daily) && prepared_climate.co2_daily
        required_co2_values = co2_daily ? local_end_day : div(local_end_day - 1, 365) + 1
        label = co2_daily ? "daily" : "annual"
        length(prepared_climate.co2) >= required_co2_values || throw(DimensionMismatch(
            "$label CO₂ forcing has $(length(prepared_climate.co2)) value(s), " *
            "but climate rows through day $local_end_day require $required_co2_values",
        ))
    elseif ndims(prepared_climate.co2) == 2
        size(prepared_climate.co2, 1) >= local_end_day || throw(DimensionMismatch(
            "daily CO₂ forcing has $(size(prepared_climate.co2, 1)) row(s), " *
            "but end_day is $local_end_day",
        ))
    else
        throw(ArgumentError("climate.co2 must be an annual vector or daily matrix"))
    end
    run_days = local_end_day - start_day + 1
    run_days <= remaining_days || throw(DimensionMismatch(
        "only $remaining_days of $(simulation.config.days) configured days remain, requested $run_days",
    ))

    if spinup && simulation.simulated_days == 0 && hasproperty(prepared_climate, :temp_spinup)
        spin_up_climbuf!(
            simulation.cft,
            prepared_climate.temp_spinup,
            simulation.climbuf;
            prec_spinup = hasproperty(prepared_climate, :prec_spinup) ?
                prepared_climate.prec_spinup : nothing,
            year_spinup = spinup_years,
        )
    end

    prescribed = isnothing(management) ?
        (prescribed_phu = nothing, prescribed_winter_type = nothing) :
        _prepare_annual_management!(simulation, management)
    _transition_range!(
        simulation, prepared_climate, start_day, local_end_day;
        reuse_output, selected_output, prescribed...,
        neural_parameters, neural_layout,
    )
    simulation.simulated_days += run_days
    return simulation
end

"""Load the `climate` object from one JLD2 file and append it to a simulation."""
function run_simulation!(
    simulation::CropSimulation,
    climate_file::AbstractString;
    kwargs...,
)
    climate = load(climate_file, "climate")
    climate isa NamedTuple || throw(ArgumentError(
        "JLD2 variable `climate` must be a NamedTuple, got $(typeof(climate))",
    ))
    return run_simulation!(simulation, climate; kwargs...)
end

"""
    run_simulation!(simulation, climate_blocks;
                    spinup=true, spinup_years=1, output_stream=nothing)

Run an ordered collection of climate `NamedTuple`s or JLD2 file paths without
concatenating them in memory. The first block may initialize the climate
buffer; subsequent blocks inherit all crop and soil state.
"""
function run_simulation!(
    simulation::CropSimulation,
    climate_blocks::AbstractVector;
    spinup::Bool = true,
    spinup_years::Integer = 1,
    output_stream::Union{Nothing, OutputStream} = nothing,
    management_blocks::Union{Nothing, AbstractVector} = nothing,
    neural_parameters = nothing,
    neural_layout = nothing,
)
    isempty(climate_blocks) && throw(ArgumentError("climate_blocks must not be empty"))
    management_blocks === nothing || length(management_blocks) == length(climate_blocks) ||
        throw(DimensionMismatch("management_blocks must match climate_blocks"))
    if output_stream !== nothing
        all(isnothing, values(simulation.diagnostics)) || throw(ArgumentError(
            "streamed global runs require initialize_simulation(...; diagnostics=false)",
        ))
        _output_timeseries_empty(simulation.output) || throw(ArgumentError(
            "streaming must start with empty output time series",
        ))
    end
    selected_output = output_stream === nothing ? nothing : Set(
        (variable.group, variable.field) for variable in output_stream.variables
    )
    try
        for (block_index, block) in pairs(climate_blocks)
            first_day = simulation.simulated_days + 1
            run_simulation!(
                simulation, block;
                spinup = spinup && block_index == firstindex(climate_blocks),
                spinup_years = spinup_years,
                reuse_output = output_stream !== nothing,
                selected_output = selected_output,
                management = management_blocks === nothing ? nothing :
                    management_blocks[block_index],
                neural_parameters,
                neural_layout,
            )
            if output_stream !== nothing
                consume_output!(
                    output_stream, simulation.output, first_day;
                    rows = simulation.simulated_days - first_day + 1,
                )
            end
        end
    finally
        applicable(close, climate_blocks) && close(climate_blocks)
    end
    if output_stream !== nothing && simulation.simulated_days == simulation.config.days
        finish_output_stream!(output_stream, simulation.simulated_days)
    end
    output_stream !== nothing && clear_output_timeseries!(simulation.output)
    return simulation
end

const _CHECKPOINT_FORMAT_VERSION = 7

_checkpoint_snapshot(values::AbstractArray) = Array(values)
_checkpoint_snapshot(values::NamedTuple) = map(_checkpoint_snapshot, values)
_checkpoint_snapshot(::Nothing) = nothing
function _checkpoint_snapshot(value)
    names = fieldnames(typeof(value))
    return NamedTuple{names}(map(
        name -> _checkpoint_snapshot(getfield(value, name)), names,
    ))
end

function _checkpoint_fields(value, names::Tuple)
    return NamedTuple{names}(map(
        name -> _checkpoint_snapshot(getproperty(value, name)), names,
    ))
end

function _restore_checkpoint_fields!(target, snapshot)
    for name in keys(snapshot)
        destination = getproperty(target, name)
        source = getproperty(snapshot, name)
        if destination isa AbstractArray
            size(destination) == size(source) || throw(DimensionMismatch(
                "checkpoint field $name has size $(size(source)); " *
                "target has size $(size(destination))",
            ))
            copyto!(destination, source)
        else
            _restore_checkpoint_fields!(destination, source)
        end
    end
    return nothing
end

function _restore_checkpoint_output!(target, snapshot, device)
    for name in keys(snapshot)
        destination = getproperty(target, name)
        source = getproperty(snapshot, name)
        if destination isa AbstractArray
            setproperty!(target, name, device(copy(source)))
        else
            _restore_checkpoint_output!(destination, source, device)
        end
    end
    return nothing
end

function _simulation_checkpoint(simulation::CropSimulation)
    cells = length(simulation.managed_land.latitude)
    return (
        format_version = _CHECKPOINT_FORMAT_VERSION,
        metadata = (
            state_schema_version = 2,
            precision = string(simulation.config.T),
            cells = cells,
            cell_ids = copy(simulation.config.execution.domain.cell_ids),
            configured_days = simulation.config.days,
            cft_id = simulation.cft.name,
            photosynthetic_pathway = simulation.cft.path,
            irrigation = simulation.config.irrigation,
            manure = simulation.config.manure,
            fertilizer = simulation.config.fertilizer,
            with_tillage = simulation.config.with_tillage,
            crop_resp_fix = simulation.config.crop_resp_fix,
            nitrogen_limit_vcmax = simulation.config.nitrogen_limit_vcmax,
            sowing_mode = simulation.config.sowing_mode,
            parameter_fingerprint = _checkpoint_fingerprint((
                cft = simulation.cft,
                model_parameters = simulation.model_parameters,
            )),
        ),
        simulated_days = simulation.simulated_days,
        cft = simulation.cft,
        model_parameters = simulation.model_parameters,
        state = (
            prognostic = _checkpoint_snapshot(simulation.state.prognostic),
            inputs = (
                crop = _checkpoint_snapshot(simulation.state.inputs.crop),
                soil = _checkpoint_snapshot(simulation.state.inputs.soil),
                management = _checkpoint_snapshot(simulation.state.inputs.management),
            ),
            output = _checkpoint_snapshot(simulation.state.output),
        ),
        diagnostics = _checkpoint_snapshot(simulation.diagnostics),
    )
end

"""Return a stable digest for immutable checkpoint metadata."""
function _checkpoint_fingerprint(value)
    return bytes2hex(sha256(codeunits(repr(value))))
end

"""
    save_checkpoint(path, simulation)

Write a backend-independent checkpoint at a completed daily boundary. Arrays
are stored on the host so the checkpoint can be restored into either a CPU or
CUDA simulation with the same precision, dimensions, and run configuration.
"""
function save_checkpoint(path::AbstractString, simulation::CropSimulation)
    checkpoint = _simulation_checkpoint(simulation)
    jldsave(path; checkpoint = checkpoint)
    return path
end

function _validate_checkpoint_target(simulation::CropSimulation, checkpoint)
    checkpoint.format_version == _CHECKPOINT_FORMAT_VERSION || throw(ArgumentError(
        "unsupported checkpoint format version $(checkpoint.format_version)",
    ))
    simulation.simulated_days == 0 || throw(ArgumentError(
        "restore_checkpoint! requires a newly initialized simulation",
    ))
    metadata = checkpoint.metadata
    checks = (
        ("state schema version", metadata.state_schema_version, 2),
        ("precision", metadata.precision, string(simulation.config.T)),
        ("cell count", metadata.cells, length(simulation.managed_land.latitude)),
        ("cell ids", metadata.cell_ids, simulation.config.execution.domain.cell_ids),
        ("configured days", metadata.configured_days, simulation.config.days),
        ("CFT id", metadata.cft_id, simulation.cft.name),
        ("photosynthetic pathway", metadata.photosynthetic_pathway, simulation.cft.path),
        ("irrigation", metadata.irrigation, simulation.config.irrigation),
        ("manure", metadata.manure, simulation.config.manure),
        ("fertilizer", metadata.fertilizer, simulation.config.fertilizer),
        ("tillage", metadata.with_tillage, simulation.config.with_tillage),
        ("fixed crop respiration ratio", metadata.crop_resp_fix,
         simulation.config.crop_resp_fix),
        ("nitrogen Vcmax limitation", metadata.nitrogen_limit_vcmax,
         simulation.config.nitrogen_limit_vcmax),
        ("sowing mode", metadata.sowing_mode, simulation.config.sowing_mode),
        ("parameter fingerprint", metadata.parameter_fingerprint,
         _checkpoint_fingerprint((
             cft = simulation.cft,
             model_parameters = simulation.model_parameters,
         ))),
    )
    for (label, saved, target) in checks
        saved == target || throw(ArgumentError(
            "checkpoint $label is $saved; target simulation uses $target",
        ))
    end
    0 <= checkpoint.simulated_days <= simulation.config.days || throw(ArgumentError(
        "checkpoint simulated_days is outside the configured simulation range",
    ))
    return nothing
end

"""
    restore_checkpoint!(simulation, path)

Restore a checkpoint into a newly initialized, configuration-compatible
simulation. Runtime arrays are copied to the target simulation's active backend.
"""
function restore_checkpoint!(simulation::CropSimulation, path::AbstractString)
    checkpoint = load(path, "checkpoint")
    _validate_checkpoint_target(simulation, checkpoint)

    simulation.processes = ProcessModules(checkpoint.cft, checkpoint.model_parameters)
    _restore_checkpoint_fields!(
        simulation.state.prognostic, checkpoint.state.prognostic,
    )
    _restore_checkpoint_fields!(simulation.state.inputs.crop, checkpoint.state.inputs.crop)
    _restore_checkpoint_fields!(simulation.state.inputs.soil, checkpoint.state.inputs.soil)
    _restore_checkpoint_fields!(
        simulation.state.inputs.management, checkpoint.state.inputs.management,
    )
    _restore_checkpoint_output!(
        simulation.state.output, checkpoint.state.output, simulation.config.device,
    )
    for name in _SIMULATION_DIAGNOSTIC_PROPERTIES
        target = getproperty(simulation.diagnostics, name)
        saved = getproperty(checkpoint.diagnostics, name)
        if target === nothing || saved === nothing
            target === saved || throw(ArgumentError(
                "checkpoint diagnostics setting does not match target simulation",
            ))
        else
            _restore_checkpoint_fields!(target, saved)
        end
    end
    simulation.simulated_days = checkpoint.simulated_days
    return simulation
end

_simulation_host(values) = Array(values)

"""Return a compact, backend-independent summary of a completed simulation."""
function simulation_summary(simulation::CropSimulation)
    simulation.simulated_days > 0 || throw(ArgumentError(
        "run the simulation before requesting its summary",
    ))
    npp = _simulation_host(simulation.output.crop.npp)
    lai = _simulation_host(simulation.output.crop.lai)
    biomass = _simulation_host(simulation.output.crop.biomass)
    water = simulation.water_balance
    nitrogen = simulation.nitrogen_balance
    carbon = simulation.carbon_balance
    thermal = simulation.thermal_balance

    balance_summary = water === nothing ? nothing : (
        maximum_absolute_residual = maximum(abs, _simulation_host(water.residual)),
        cumulative_residual = sum(_simulation_host(water.residual)),
    )
    nitrogen_summary = nitrogen === nothing ? nothing : (
        maximum_absolute_residual = maximum(abs, _simulation_host(nitrogen.residual)),
        cumulative_residual = sum(_simulation_host(nitrogen.residual)),
        cumulative_leaching_loss = sum(_simulation_host(nitrogen.leaching_loss)),
    )
    carbon_summary = carbon === nothing ? nothing : (
        maximum_absolute_residual = maximum(abs, _simulation_host(carbon.residual)),
        cumulative_residual = sum(_simulation_host(carbon.residual)),
        cumulative_npp = sum(_simulation_host(carbon.net_primary_production)),
    )
    thermal_summary = thermal === nothing ? nothing : (
        maximum_absolute_energy_residual = maximum(abs, _simulation_host(thermal.energy_residual)),
        maximum_absolute_percolation_energy_residual = maximum(
            abs, _simulation_host(thermal.percolation_energy_residual),
        ),
    )
    return (
        precision = simulation.config.T,
        cells = size(npp, 2),
        simulated_days = simulation.simulated_days,
        crop = (
            cumulative_npp = sum(npp),
            maximum_lai = maximum(lai),
            maximum_biomass = maximum(biomass),
        ),
        water = balance_summary,
        nitrogen = nitrogen_summary,
        carbon = carbon_summary,
        thermal = thermal_summary,
    )
end
