using TOML
using Dates
using SHA

include(joinpath(@__DIR__, "run_global_wheat_cpu.jl"))

function requested_crop_systems(config)
    cfts = get(config, "cfts", Dict{String, Any}())
    requested = get(cfts, "cft_ids", [1])
    cft_ids = requested isa AbstractString && lowercase(requested) == "all" ?
        collect(1:length(CFTS)) : Int.(requested)
    isempty(cft_ids) && error("cfts.cft_ids must not be empty")
    all(id -> 1 <= id <= length(CFTS), cft_ids) ||
        error("cfts.cft_ids must be in 1:$(length(CFTS))")
    allunique(cft_ids) || error("cfts.cft_ids must be unique")
    systems = Symbol.(lowercase.(String.(get(cfts, "water_systems", ["rainfed"]))))
    all(system -> system in (:rainfed, :irrigated), systems) || error(
        "cfts.water_systems must contain rainfed and/or irrigated",
    )
    allunique(systems) || error("cfts.water_systems must be unique")
    return [(cft_id, system === :irrigated) for cft_id in cft_ids for system in systems]
end

batch_name(cft_id, irrigated) = "cft_$(lpad(cft_id, 2, '0'))_" *
    (irrigated ? "irrigated" : "rainfed")

const _FREE_WARMUP_RUN_KEYS = (
    "warmup_minimum_years",
    "warmup_maximum_years",
    "warmup_consecutive_years",
    "warmup_relative_tolerance",
    "warmup_pool_fraction_tolerance",
    "warmup_required_converged_fraction",
    "warmup_cache_climate",
    "require_warmup_convergence",
)

const _CALIBRATION_RUN_OPTIONS = Dict{String, Any}(
    "warmup_minimum_years" => 600,
    "warmup_maximum_years" => 600,
    "warmup_target_constrained" => true,
    # Calibration always writes its allocation product. The free warm-up stage
    # is responsible for deciding whether the resulting state is production-ready.
    "require_warmup_convergence" => false,
)

const _CFT_BATCH_MANIFEST_SCHEMA_VERSION = "1"

_file_fingerprint(path::AbstractString) = bytes2hex(sha256(read(path)))

function _repository_commit()
    try
        return readchomp(`git rev-parse --verify HEAD`)
    catch
        return "unknown"
    end
end

function _write_batch_status(path; status, name, cft_id, irrigated,
    calibration_config, production_config, production_output,
    calibration_output_directory, production_output_directory,
    allocation_path, production_config_fingerprint, kwargs...)
    payload = Dict{String, Any}(
        "schema_version" => _CFT_BATCH_MANIFEST_SCHEMA_VERSION,
        "status" => String(status),
        "name" => name,
        "cft_id" => cft_id,
        "irrigated" => irrigated,
        "calibration_config" => abspath(calibration_config),
        "production_config" => abspath(production_config),
        "production_output" => abspath(production_output),
        "calibration_output_directory" => abspath(calibration_output_directory),
        "production_output_directory" => abspath(production_output_directory),
        "pool_allocation" => abspath(allocation_path),
        "production_config_fingerprint" => production_config_fingerprint,
        "updated_at" => string(now()),
    )
    for (key, value) in pairs(kwargs)
        payload[String(key)] = value
    end
    write_report(path, payload)
    return path
end

function _resumable_batch(status_path, production_config_fingerprint)
    isfile(status_path) || return nothing
    status = TOML.parsefile(status_path)
    get(status, "schema_version", "") == _CFT_BATCH_MANIFEST_SCHEMA_VERSION || return nothing
    get(status, "status", "") == "completed" || return nothing
    get(status, "production_config_fingerprint", "") == production_config_fingerprint || return nothing
    production_output = get(status, "production_output", nothing)
    allocation_path = get(status, "pool_allocation", nothing)
    (production_output isa AbstractString && isfile(production_output)) || return nothing
    (allocation_path isa AbstractString && isfile(allocation_path)) || return nothing
    return status
end

function _resumable_calibration(status_path, production_config_fingerprint)
    isfile(status_path) || return nothing
    status = TOML.parsefile(status_path)
    get(status, "schema_version", "") == _CFT_BATCH_MANIFEST_SCHEMA_VERSION || return nothing
    status["status"] in ("calibration_completed", "failed") || return nothing
    get(status, "calibration_completed", false) || return nothing
    get(status, "production_config_fingerprint", "") == production_config_fingerprint || return nothing
    allocation_path = get(status, "pool_allocation", nothing)
    allocation_fingerprint = get(status, "allocation_fingerprint", nothing)
    (allocation_path isa AbstractString && isfile(allocation_path)) || return nothing
    allocation_fingerprint == _file_fingerprint(allocation_path) || return nothing
    return status
end

function write_batch_config(path, config; output_directory, pool_allocation = nothing,
    production::Bool, target_constrained::Bool, warmup_options = nothing)
    batch = deepcopy(config)
    paths = batch["paths"]
    run = batch["run"]
    paths["output_directory"] = output_directory
    if isnothing(pool_allocation)
        pop!(paths, "pool_allocation", nothing)
    else
        paths["pool_allocation"] = pool_allocation
    end
    run["production"] = production
    run["calibration_only"] = !production
    run["warmup_target_constrained"] = target_constrained
    if !production
        merge!(run, _CALIBRATION_RUN_OPTIONS)
    elseif !isnothing(warmup_options)
        for key in _FREE_WARMUP_RUN_KEYS
            haskey(warmup_options, key) && (run[key] = warmup_options[key])
        end
    end
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, batch; sorted = true)
    end
    return path
end

function run_global_cfts(
    config_path;
    backend_override = :cpu,
    cft_id::Union{Nothing, Integer} = nothing,
    irrigated::Union{Nothing, Bool} = nothing,
    output_directory_override::Union{Nothing, AbstractString} = nothing,
    partition_rank::Integer = 0,
    partition_count::Integer = 1,
    warmup_convergence_reducer = nothing,
    resume_consensus = nothing,
)
    config = TOML.parsefile(config_path)
    haskey(config["paths"], "pool_allocation") && error(
        "the CFT workflow derives one CFT- and water-system-specific allocation per batch; " *
        "do not set paths.pool_allocation in the shared configuration",
    )
    systems = requested_crop_systems(config)
    if !isnothing(cft_id)
        1 <= cft_id <= length(CFTS) || error("cft_id must be in 1:$(length(CFTS))")
        isnothing(irrigated) && error("a single CFT run must specify rainfed or irrigated")
        systems = [(Int(cft_id), irrigated)]
    end
    output_directory = isnothing(output_directory_override) ?
        abspath(config["paths"]["output_directory"]) : abspath(output_directory_override)
    haskey(config, "free_warmup") && haskey(config, "allocation_validation") && error(
        "use [free_warmup], not both [free_warmup] and legacy [allocation_validation]",
    )
    free_warmup_options = get(config, "free_warmup",
        get(config, "allocation_validation", nothing))
    calibration_only = Bool(get(config["run"], "calibration_only", false))
    simulation_years = configured_simulation_years(config)
    resume_completed_batches = Bool(get(config["run"], "resume_completed_batches", true))
    batch_manifest = Dict{String, Any}[]
    for (cft_id, irrigated) in systems
        name = batch_name(cft_id, irrigated)
        println("running $name")
        batch_directory = joinpath(output_directory, "batches", name)
        calibration_directory = joinpath(batch_directory, "calibration")
        calibration_config = write_batch_config(
            joinpath(batch_directory, "calibration.toml"), config;
            output_directory = calibration_directory,
            production = false,
            target_constrained = true,
        )
        production_directory = joinpath(batch_directory, "production")
        allocation_path = joinpath(calibration_directory, "warmup_soil_pool_allocation.nc")
        production_config = write_batch_config(
            joinpath(batch_directory, "production.toml"), config;
            output_directory = production_directory,
            pool_allocation = allocation_path,
            production = true,
            target_constrained = false,
            warmup_options = free_warmup_options,
        )
        production_output = joinpath(
            production_directory,
            "global_wheat_$(first(simulation_years))_$(last(simulation_years)).nc",
        )
        production_config_fingerprint = _file_fingerprint(production_config)
        status_path = joinpath(batch_directory, "batch_status.toml")
        resumed = resume_completed_batches ?
            _resumable_batch(status_path, production_config_fingerprint) : nothing
        if resume_completed_batches && resume_consensus !== nothing
            globally_completed = resume_consensus(resumed !== nothing)
            globally_completed isa Bool || error("resume consensus must return Bool")
            resumed = globally_completed ? resumed : nothing
        end
        if resumed !== nothing
            println("resuming completed $name")
            push!(batch_manifest, Dict{String, Any}(String(k) => v for (k, v) in pairs(resumed)))
            continue
        end

        calibration_resume = resume_completed_batches ?
            _resumable_calibration(status_path, production_config_fingerprint) : nothing
        if resume_completed_batches && resume_consensus !== nothing
            globally_calibrated = resume_consensus(calibration_resume !== nothing)
            globally_calibrated isa Bool || error("resume consensus must return Bool")
            calibration_resume = globally_calibrated ? calibration_resume : nothing
        end
        calibration_resume === nothing && _write_batch_status(
            status_path;
            status = "running", name, cft_id, irrigated,
            calibration_config, production_config, production_output,
            calibration_output_directory = calibration_directory,
            production_output_directory = production_directory,
            allocation_path, production_config_fingerprint,
            partition_rank = Int(partition_rank),
            partition_count = Int(partition_count),
        )
        try
            calibration_output_directory = calibration_directory
            calibration_completed = false
            if calibration_resume === nothing
                calibration = run_global_wheat(
                    calibration_config;
                    backend_override,
                    cft_id,
                    irrigated,
                    partition_rank,
                    partition_count,
                    warmup_convergence_reducer,
                )
                isfile(allocation_path) || error("$name did not write a soil-pool allocation")
                allocation_fingerprint = _file_fingerprint(allocation_path)
                calibration_output_directory = calibration.output_directory
                _write_batch_status(
                    status_path;
                    status = "calibration_completed", name, cft_id, irrigated,
                    calibration_config, production_config, production_output,
                    calibration_output_directory,
                    production_output_directory = production_directory,
                    allocation_path, production_config_fingerprint,
                    allocation_fingerprint,
                    calibration_completed = true,
                    partition_rank = Int(partition_rank),
                    partition_count = Int(partition_count),
                )
                calibration_completed = true
            else
                println("resuming production after completed calibration for $name")
                allocation_fingerprint = calibration_resume["allocation_fingerprint"]
                calibration_output_directory = calibration_resume["calibration_output_directory"]
                calibration_completed = true
            end
            if calibration_only
                calibration_status = TOML.parsefile(status_path)
                push!(batch_manifest, Dict{String, Any}(
                    String(key) => value for (key, value) in pairs(calibration_status)
                ))
                continue
            end
            result = run_global_wheat(
                production_config;
                backend_override,
                cft_id,
                irrigated,
                partition_rank,
                partition_count,
                warmup_convergence_reducer,
            )
            isfile(production_output) || error("$name did not write production output")
            status = Dict{String, Any}(
                "schema_version" => _CFT_BATCH_MANIFEST_SCHEMA_VERSION,
                "status" => "completed",
                "name" => name,
                "cft_id" => cft_id,
                "irrigated" => irrigated,
                "calibration_config" => abspath(calibration_config),
                "production_config" => abspath(production_config),
                "production_output" => abspath(production_output),
                "calibration_output_directory" => abspath(calibration_output_directory),
                "production_output_directory" => abspath(result.output_directory),
                "pool_allocation" => abspath(allocation_path),
                "production_config_fingerprint" => production_config_fingerprint,
                "allocation_fingerprint" => allocation_fingerprint,
                "production_output_bytes" => filesize(production_output),
                "production_warmup_years" => result.warmup_years,
                "production_warmup_converged" => result.warmup_converged,
                "partition_rank" => Int(partition_rank),
                "partition_count" => Int(partition_count),
                "updated_at" => string(now()),
            )
            write_report(status_path, status)
            push!(batch_manifest, status)
        catch exception
            _write_batch_status(
                status_path;
                status = "failed", name, cft_id, irrigated,
                calibration_config, production_config, production_output,
                calibration_output_directory = calibration_directory,
                production_output_directory = production_directory,
                allocation_path, production_config_fingerprint,
                calibration_completed,
                partition_rank = Int(partition_rank),
                partition_count = Int(partition_count),
                error = sprint(showerror, exception),
            )
            rethrow()
        end
    end
    manifest_path = joinpath(output_directory, "cft_batch_manifest.toml")
    write_report(manifest_path, Dict(
        "schema_version" => _CFT_BATCH_MANIFEST_SCHEMA_VERSION,
        "repository_commit" => _repository_commit(),
        "config_path" => abspath(config_path),
        "config_fingerprint" => _file_fingerprint(config_path),
        "simulation_start_year" => first(simulation_years),
        "simulation_end_year" => last(simulation_years),
        "partition_rank" => Int(partition_rank),
        "partition_count" => Int(partition_count),
        "created_at" => string(now()),
        "batches" => batch_manifest,
    ))
    return manifest_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (1, 3) || error(
        "usage: run_global_cfts_cpu.jl CONFIG_TOML [CFT_ID rainfed|irrigated]",
    )
    cft_id = length(ARGS) == 3 ? parse(Int, ARGS[2]) : nothing
    irrigated = length(ARGS) == 3 ?
        Symbol(lowercase(ARGS[3])) === :irrigated : nothing
    length(ARGS) == 3 && Symbol(lowercase(ARGS[3])) in (:rainfed, :irrigated) ||
        length(ARGS) == 1 || error("water system must be rainfed or irrigated")
    println(run_global_cfts(
        abspath(ARGS[1]); backend_override = :cpu, cft_id, irrigated,
    ))
end
