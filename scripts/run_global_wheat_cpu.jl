using NeuralCrop
using Dates
using NCDatasets
using TOML

if !isdefined(@__MODULE__, :NeuralCropData)
    include(joinpath(@__DIR__, "..", "lib", "NeuralCropData", "src", "NeuralCropData.jl"))
end
import .NeuralCropData: CFTRegistry, DATA_SCHEMA_VERSION, DatasetCatalog, DatasetSpec,
    HWSD_CN_PREPROCESSING_VERSION, PatchDomain, SoilCNTargets, SoilPoolAllocation,
    build_crop_mask, cft_index, climate_blocks, climate_days, climate_forcings,
    combine_patch_domains, dataset, expand_to_grid, management_schedule,
    read_climate_block, read_grid, read_management, read_soil_cn_targets,
    read_soil_data, read_soil_pool_allocation, select_cells, write_soil_pool_allocation

function execution_backend(config; override = nothing)
    configured = Symbol(lowercase(String(get(config["run"], "backend", "cpu"))))
    backend = isnothing(override) ? configured : Symbol(override)
    backend in (:cpu, :cuda) || error("run.backend must be cpu or cuda")
    backend === :cpu && return (name = :cpu, device = identity, device_id = nothing)

    CUDA.functional() || error("CUDA backend requested, but CUDA.functional() is false")
    device_id = Int(get(config["run"], "device_id", 0))
    CUDA.device!(device_id)
    CUDA.allowscalar(false)
    return (name = :cuda, device = CUDA.CuArray, device_id)
end

function catalog_from_config(config)
    paths = config["paths"]
    input_directory = abspath(paths["input_directory"])
    management_directory = joinpath(input_directory, "management")
    climate_directory = joinpath(input_directory, "climate")
    soil_directory = joinpath(input_directory, "soil")
    climate = get(config, "climate", Dict{String, Any}())
    climate_path(name, default) = begin
        path = String(get(climate, name, default))
        isabspath(path) ? path : joinpath(climate_directory, path)
    end
    management_config = get(config, "management", Dict{String, Any}())
    management_path(name, default) = begin
        path = String(get(management_config, name, default))
        isabspath(path) ? path : joinpath(management_directory, path)
    end
    registry = CFTRegistry(collect(1:12), collect(CFT_NAMES))
    bands = (rainfed = Int32.(1:12), irrigated = Int32.(13:24))
    multi_cft(path, variable; residue = false) = DatasetSpec(
        path, variable;
        rainfed_bands = bands.rainfed,
        irrigated_bands = residue ? bands.rainfed : bands.irrigated,
    )
    management = Dict{Symbol, DatasetSpec}(
        :landuse => multi_cft(
            management_path("landuse_file", "landuse_24cfts_2015.nc"), "landfrac",
        ),
        :sowing_date => multi_cft(
            management_path("sowing_date_file", "sdate_24cfts_2015.nc"), "sdate",
        ),
        :phu => multi_cft(
            management_path("phu_file", "phu_24cfts_2015.nc"), "phusum",
        ),
        :fertilizer => multi_cft(
            management_path("fertilizer_file", "fertilizer_24cfts_2015.nc"), "fertilizer",
        ),
        :manure => multi_cft(
            management_path("manure_file", "manure_24cfts_2015.nc"), "manure",
        ),
        :residue_fraction => multi_cft(
            management_path("residue_fraction_file", "residue_12cfts_2015.nc"),
            "residuefrac"; residue = true,
        ),
    )
    has_no3_deposition = haskey(climate, "no3_deposition_file")
    has_nh4_deposition = haskey(climate, "nh4_deposition_file")
    has_no3_deposition == has_nh4_deposition || error(
        "configure both climate.no3_deposition_file and climate.nh4_deposition_file",
    )
    climate_specs = Dict{Symbol, DatasetSpec}(
        :grid => DatasetSpec(joinpath(soil_directory, "grid.nc"), "cellid"),
        :soilcode => DatasetSpec(joinpath(soil_directory, "soil_30arcmin_13_types.nc"), "soilcode"),
        :soilph => DatasetSpec(joinpath(soil_directory, "soil_pH30arcmin.nc"), "soilph"),
        :temp => DatasetSpec(climate_path("temperature_file", "temp_2015_2016.nc"), "temp"),
        :prec => DatasetSpec(climate_path("precipitation_file", "prec_2015_2016.nc"), "prec"),
        :lwnet => DatasetSpec(climate_path("longwave_file", "lwnet_2015_2016.nc"), "lwnet"),
        :swdown => DatasetSpec(climate_path("shortwave_file", "swdown_2015_2016.nc"), "swdown"),
        :co2 => DatasetSpec(climate_path("co2_file", "co2_2015_2016.txt"), "co2"),
    )
    if has_no3_deposition
        climate_specs[:no3_deposition] = DatasetSpec(
            climate_path("no3_deposition_file", ""),
            String(get(climate, "no3_deposition_variable", "noy"));
            units = "g/m2/day",
        )
        climate_specs[:nh4_deposition] = DatasetSpec(
            climate_path("nh4_deposition_file", ""),
            String(get(climate, "nh4_deposition_variable", "nhx"));
            units = "g/m2/day",
        )
    end
    return DatasetCatalog(
        merge(management, climate_specs),
        registry,
    )
end

function selected_targets(targets, selection)
    position = Dict(cell_id => index for (index, cell_id) in pairs(targets.selection.cell_ids))
    indices = [get(position, cell_id, 0) for cell_id in selection.cell_ids]
    all(>(0), indices) || error("HWSD targets do not cover every selected wheat cell")
    return SoilCNTargets(
        selection,
        targets.layer_bounds,
        targets.soil_organic_carbon[:, indices],
        targets.total_nitrogen[:, indices],
        targets.coverage[:, indices],
        targets.uncertain[:, indices],
        targets.provenance,
    )
end

function selected_pool_allocation(allocation, selection)
    allocation === nothing && return nothing
    position = Dict(cell_id => index for (index, cell_id) in pairs(allocation.selection.cell_ids))
    indices = [get(position, cell_id, 0) for cell_id in selection.cell_ids]
    all(>(0), indices) || error("pool allocation does not cover every selected wheat cell")
    return SoilPoolAllocation(
        selection,
        allocation.fast_carbon_fraction[:, indices],
        allocation.fast_nitrogen_fraction[:, indices],
        allocation.c_shift_fast[:, indices],
        allocation.c_shift_slow[:, indices];
        cft_id = allocation.cft_id,
        irrigated = allocation.irrigated,
        provenance = allocation.provenance,
    )
end

function load_hwsd_targets(paths, selection)
    if haskey(paths, "hwsd_targets") && isfile(abspath(paths["hwsd_targets"]))
        return read_soil_cn_targets(abspath(paths["hwsd_targets"]); T = Float32)
    end
    haskey(paths, "hwsd_profile_directory") || error(
        "provide an existing hwsd_targets file or hwsd_profile_directory",
    )
    directory = abspath(paths["hwsd_profile_directory"])
    profiles = [
        read_soil_cn_targets(joinpath(directory, "cell_$(cell_id).nc"); T = Float32)
        for cell_id in selection.cell_ids
    ]
    return SoilCNTargets(
        selection,
        first(profiles).layer_bounds,
        reduce(hcat, (profile.soil_organic_carbon for profile in profiles)),
        reduce(hcat, (profile.total_nitrogen for profile in profiles)),
        reduce(hcat, (profile.coverage for profile in profiles)),
        BitMatrix(reduce(hcat, (profile.uncertain for profile in profiles))),
        (
            schema_version = DATA_SCHEMA_VERSION,
            preprocessing_version = HWSD_CN_PREPROCESSING_VERSION,
            source_version = "HWSD v2.01 per-cell profiles",
            deep_rule = :extend_deepest_density,
        ),
    )
end

function management_source_years(config, simulation_years)
    management = config["management"]
    mode = Symbol(lowercase(String(get(management, "mode", "fixed"))))
    mode in (:fixed, :transient) || error("management mode must be fixed or transient")
    if mode === :fixed
        fixed_year = Int(get(management, "fixed_year", 2015))
        return fill(fixed_year, length(simulation_years))
    end
    return collect(Int, simulation_years)
end

"""Resolve an optional fixed year for atmospheric nitrogen deposition."""
function nitrogen_deposition_source_year(config)
    climate = get(config, "climate", Dict{String, Any}())
    management = config["management"]
    default_mode = get(management, "mode", "fixed")
    mode = Symbol(lowercase(String(get(
        climate, "nitrogen_deposition_mode", default_mode,
    ))))
    mode in (:fixed, :transient) || error(
        "climate.nitrogen_deposition_mode must be fixed or transient",
    )
    mode === :transient && return nothing
    return Int(get(
        climate,
        "nitrogen_deposition_fixed_year",
        get(management, "fixed_year", 2015),
    ))
end

"""Resolve production years, or the fixed-management anchor for calibration only."""
function configured_simulation_years(config)
    run = config["run"]
    has_start = haskey(run, "simulation_start_year")
    has_end = haskey(run, "simulation_end_year")
    has_start == has_end || error(
        "provide both simulation_start_year and simulation_end_year",
    )
    if has_start
        start_year = Int(run["simulation_start_year"])
        end_year = Int(run["simulation_end_year"])
        start_year <= end_year || error(
            "simulation_start_year must not exceed simulation_end_year",
        )
        return start_year:end_year
    end

    calibration_only = Bool(get(run, "calibration_only", false)) ||
        !Bool(get(run, "production", true))
    calibration_only || error(
        "production runs require simulation_start_year and simulation_end_year",
    )
    management = config["management"]
    mode = Symbol(lowercase(String(get(management, "mode", "fixed"))))
    mode === :fixed || error(
        "calibration without simulation years requires management.mode = \"fixed\"",
    )
    fixed_year = Int(get(management, "fixed_year", 2015))
    return fixed_year:fixed_year
end

function read_management_schedule(
    catalog, grid, selection, active, source_years, simulation_years, config, cft_id;
    irrigated::Bool,
)
    sowing_date = read_management(
        catalog, :sowing_date, grid, cft_id;
        selection, active, simulation_years = source_years, T = Float32, irrigated,
    )
    phu = read_management(
        catalog, :phu, grid, cft_id;
        selection, active, simulation_years = source_years, T = Float32, irrigated,
    )
    fertilizer = read_management(
        catalog, :fertilizer, grid, cft_id;
        selection, simulation_years = source_years, T = Float32, irrigated,
    )
    manure = read_management(
        catalog, :manure, grid, cft_id;
        selection, simulation_years = source_years, T = Float32, irrigated,
    )
    residue = read_management(
        catalog, :residue_fraction, grid, cft_id;
        selection, simulation_years = source_years, T = Float32, irrigated,
    )
    management = config["management"]
    schedule = management_schedule(
        ; sowing_date, phu, fertilizer, manure, residue_fraction = residue,
        fertilizer_mode = Symbol(management["fertilizer"]),
        manure_enabled = management["manure"],
    )
    return merge(schedule, (years = Int32.(simulation_years),))
end

annual_management(schedule, index) = (
    sdate = vec(schedule.sdate[index, :]),
    phu = vec(schedule.phu[index, :]),
    manure = vec(schedule.manure[index, :]),
    fertilizer = vec(schedule.fertilizer[index, :]),
    residuefrac = vec(schedule.residuefrac[index, :]),
)

function model_inputs(
    grid, selection, hwsd_targets, catalog, management;
    pool_allocation = nothing,
)
    crop = annual_management(management, 1)
    soil = read_soil_data(catalog, grid; selection)
    targets = selected_targets(hwsd_targets, selection)
    allocation = selected_pool_allocation(pool_allocation, selection)
    initial_state = soil_initial_state(targets, soil; allocation)
    return model_initial_data(grid, soil, crop, initial_state)
end

function create_simulation(initial_data, selection, config, days, device, cft_id;
    irrigated::Bool, diagnostics)
    management = config["management"]
    sowing_mode = Symbol(get(management, "sowing_mode", "prescribed_sdate"))
    return initialize_simulation(
        crop_cft(cft_id), initial_data;
        days,
        indices = collect(1:length(selection.cell_ids)),
        cell_ids = selection.cell_ids,
        device,
        T = Float32,
        diagnostics,
        irrigation = irrigated,
        manure = management["manure"],
        fertilizer = Symbol(management["fertilizer"]),
        with_tillage = management["with_tillage"],
        crop_resp_fix = Bool(get(config["run"], "crop_resp_fix", false)),
        nitrogen_limit_vcmax = false,
        sowing_mode,
        # LPJmL keeps V_req fixed once prescribed crop dates/PHU are fixed.
        freeze_vernalization_requirement = Symbol(management["mode"]) === :fixed &&
            sowing_mode === :prescribed_sdate,
    )
end

function state_equal(left, right)
    left isa AbstractArray && return Array(left) == Array(right)
    typeof(left) == typeof(right) || return false
    fields = fieldnames(typeof(left))
    isempty(fields) && return isequal(left, right)
    return all(name -> state_equal(getfield(left, name), getfield(right, name)), fields)
end

function run_output_block!(simulation, forcing, management, stream)
    first_day = simulation.simulated_days + 1
    selected = Set((variable.group, variable.field) for variable in stream.variables)
    run_simulation!(
        simulation, forcing;
        spinup = false, reuse_output = true, selected_output = selected, management,
    )
    rows = simulation.simulated_days - first_day + 1
    consume_output!(stream, simulation.output, first_day; rows)
    clear_output_timeseries!(simulation.output)
    return simulation
end

production_output_variables() = [
    OutputVariable(:crop, :gpp; reduction = :sum),
    OutputVariable(:crop, :npp; reduction = :sum),
    OutputVariable(:soil, :ecosystem_respiration; reduction = :sum),
    OutputVariable(:soil, :heterotrophic_respiration; reduction = :sum),
    OutputVariable(:soil, :evapotranspiration; reduction = :sum),
    OutputVariable(:crop, :lai; reduction = :mean),
    OutputVariable(:crop, :yield),
]

function write_reconstructed_output(path, grid, selection, chunks, years)
    length(chunks) == length(years) || error(
        "expected $(length(years)) annual output chunks, got $(length(chunks))",
    )
    names = sort!(collect(keys(chunks[1].values)); by = string)
    isfile(path) && rm(path; force = true)
    NCDataset(path, "c") do dataset
        defDim(dataset, "longitude", length(grid.longitude))
        defDim(dataset, "latitude", length(grid.latitude))
        defDim(dataset, "time", length(years))
        defVar(dataset, "longitude", Float32, ("longitude",))[:] = grid.longitude
        defVar(dataset, "latitude", Float32, ("latitude",))[:] = grid.latitude
        defVar(dataset, "time", Int32, ("time",))[:] = Int32.(years)
        defVar(dataset, "cellid", Int32, ("longitude", "latitude"))[:, :] = grid.cellid
        for name in names
            values = Array{Float32}(
                undef, length(grid.longitude), length(grid.latitude), length(years),
            )
            for year_index in eachindex(years)
                values[:, :, year_index] = expand_to_grid(
                    vec(chunks[year_index].values[name]), grid; selection,
                )
            end
            variable = defVar(dataset, String(name), Float32, ("longitude", "latitude", "time"))
            variable[:, :, :] = values
        end
    end
    return path
end

toml_value(value::NamedTuple) = Dict(string(name) => toml_value(item) for (name, item) in pairs(value))
toml_value(value::AbstractDict) = Dict(string(name) => toml_value(item) for (name, item) in pairs(value))
toml_value(value::AbstractArray) = toml_value.(collect(value))
toml_value(value::DataType) = string(value)
toml_value(value::Symbol) = string(value)
toml_value(value) = value

function write_report(path, report)
    open(path, "w") do io
        TOML.print(io, toml_value(report); sorted = true)
    end
    return path
end

"""Write per-cell same-phase warm-up diagnostics for offline convergence QC."""
function write_warmup_cell_audit(path, grid, selection, warmup)
    cells = length(selection.cell_ids)
    length(warmup.consecutive_stable_years) == cells || error(
        "warm-up convergence audit does not match the selected cells",
    )
    comparison_year = max(0, warmup.years - warmup.forcing_years)
    previous(name) = comparison_year == 0 ?
        getproperty(warmup.initial_soil, name) :
        vec(getproperty(warmup.soil, name)[comparison_year, :])
    current(name) = vec(getproperty(warmup.soil, name)[end, :])
    relative_change(name) = (current(name) .- previous(name)) ./
        max.(abs.(previous(name)), eps(Float32))
    fast_fraction(fast, slow) = fast ./ (fast .+ slow)
    current_fast_fraction = fast_fraction(current(:fast_carbon), current(:slow_carbon))
    previous_fast_fraction = fast_fraction(previous(:fast_carbon), previous(:slow_carbon))
    positions = CartesianIndices(grid.cellid)[selection.compact_indices]
    longitude = Float32[grid.longitude[position[1]] for position in positions]
    latitude = Float32[grid.latitude[position[2]] for position in positions]
    correction = warmup.target_correction
    isfile(path) && rm(path; force = true)
    NCDataset(path, "c") do dataset
        defDim(dataset, "cell", cells)
        dataset.attrib["long_name"] = "per-cell agricultural warm-up convergence audit"
        dataset.attrib["forcing_years"] = warmup.forcing_years
        dataset.attrib["warmup_years"] = warmup.years
        dataset.attrib["target_constrained"] = string(warmup.target_constrained)
        dataset.attrib["consecutive_years_required"] = warmup.consecutive_years
        defVar(dataset, "cell_id", Int32, ("cell",))[:] = selection.cell_ids
        defVar(dataset, "longitude", Float32, ("cell",))[:] = longitude
        defVar(dataset, "latitude", Float32, ("cell",))[:] = latitude
        defVar(dataset, "strict_converged", Int8, ("cell",))[:] = Int8.(
            warmup.consecutive_stable_years .>= warmup.consecutive_years,
        )
        defVar(dataset, "consecutive_stable_years", Int32, ("cell",))[:] =
            warmup.consecutive_stable_years
        defVar(dataset, "same_phase_carbon_relative_change", Float32, ("cell",))[:] =
            Float32.(relative_change(:total_carbon))
        defVar(dataset, "same_phase_nitrogen_relative_change", Float32, ("cell",))[:] =
            Float32.(relative_change(:total_nitrogen))
        defVar(dataset, "same_phase_fast_carbon_fraction_change", Float32, ("cell",))[:] =
            Float32.(current_fast_fraction .- previous_fast_fraction)
        defVar(dataset, "final_target_carbon_correction", Float32, ("cell",))[:] =
            Float32.(vec(correction.carbon[end, :]))
        defVar(dataset, "final_target_nitrogen_correction", Float32, ("cell",))[:] =
            Float32.(vec(correction.nitrogen[end, :]))
    end
    return path
end

"""Return the scale-aware closure summary for routed water enthalpy."""
function percolation_energy_closure(residual, rain, snowmelt, lateral_runoff, bottom_drainage;
                                    absolute_tolerance, relative_tolerance)
    boundary_scale = max.(
        abs.(rain) .+ abs.(snowmelt) .+ abs.(lateral_runoff) .+ abs.(bottom_drainage),
        1.0f0,
    )
    tolerance = absolute_tolerance .+ relative_tolerance .* boundary_scale
    return (
        maximum_absolute_residual = maximum(abs, residual),
        maximum_relative_residual = maximum(abs.(residual) ./ boundary_scale),
        closes = all(abs.(residual) .<= tolerance),
    )
end

function balance_report(simulation, validation)
    water = maximum(abs, Array(simulation.water_balance.residual))
    carbon = maximum(abs, Array(simulation.carbon_balance.relative_residual))
    nitrogen = maximum(abs, Array(simulation.nitrogen_balance.relative_residual))
    energy = Array(simulation.thermal_balance.energy_residual)
    column_energy = Array(simulation.thermal_balance.column_energy)
    thermal = maximum(abs.(energy) ./ max.(abs.(column_energy), 1.0f0))
    thermal_balance = simulation.thermal_balance
    percolation = percolation_energy_closure(
        Array(thermal_balance.percolation_energy_residual),
        Array(thermal_balance.rain_energy_input),
        Array(thermal_balance.snowmelt_energy_input),
        Array(thermal_balance.lateral_runoff_energy_output),
        Array(thermal_balance.bottom_drainage_energy_output);
        absolute_tolerance = validation["maximum_percolation_energy_residual"],
        relative_tolerance = get(validation, "maximum_percolation_energy_relative_residual", 5.0e-6),
    )
    water <= validation["maximum_water_residual"] || error("sampled water closure failed")
    carbon <= validation["maximum_carbon_relative_residual"] || error("sampled carbon closure failed")
    nitrogen <= validation["maximum_nitrogen_relative_residual"] || error("sampled nitrogen closure failed")
    thermal <= validation["maximum_thermal_relative_residual"] || error("sampled thermal closure failed")
    percolation.closes || error(
        "sampled percolation-energy closure failed: " *
        "maximum_absolute_residual=$(percolation.maximum_absolute_residual), " *
        "maximum_relative_residual=$(percolation.maximum_relative_residual)",
    )
    return Dict{String, Any}(
        "maximum_water_residual" => water,
        "maximum_carbon_relative_residual" => carbon,
        "maximum_nitrogen_relative_residual" => nitrogen,
        "maximum_thermal_relative_residual" => thermal,
        "maximum_percolation_energy_residual" => percolation.maximum_absolute_residual,
        "maximum_percolation_energy_relative_residual" => percolation.maximum_relative_residual,
    )
end

"""Return one balanced, contiguous, zero-based partition of a cell selection."""
function partition_cell_selection(
    selection::NeuralCropData.CellSelection,
    rank::Integer,
    count::Integer,
)
    count > 0 || throw(ArgumentError("partition count must be positive"))
    0 <= rank < count || throw(ArgumentError("partition rank must be in 0:$(count - 1)"))
    cells = length(selection.cell_ids)
    cells >= count || throw(ArgumentError(
        "cannot divide $cells selected cells among $count non-empty partitions",
    ))
    base, remainder = divrem(cells, count)
    length_for_rank = base + (rank < remainder)
    first_index = rank * base + min(rank, remainder) + 1
    indices = first_index:(first_index + length_for_rank - 1)
    return NeuralCropData.CellSelection(
        selection.compact_indices[indices],
        selection.cell_ids[indices],
    )
end

function run_global_wheat(
    config_path;
    backend_override = nothing,
    cft_id::Integer = 1,
    irrigated::Bool = false,
    output_subdirectory::Union{Nothing, AbstractString} = nothing,
    partition_rank::Integer = 0,
    partition_count::Integer = 1,
    warmup_convergence_reducer = nothing,
)
    config = TOML.parsefile(config_path)
    paths = config["paths"]
    run = config["run"]
    production_enabled = Bool(get(run, "production", true)) &&
        !Bool(get(run, "calibration_only", false))
    backend = execution_backend(config; override = backend_override)
    println(
        "execution backend: $(backend.name)" *
        (isnothing(backend.device_id) ? "" : ", device=$(backend.device_id)"),
    )
    output_directory = abspath(paths["output_directory"])
    !isnothing(output_subdirectory) && (output_directory = joinpath(output_directory, output_subdirectory))
    mkpath(output_directory)
    catalog = catalog_from_config(config)
    grid = read_grid(dataset(catalog, :grid); T = Float32)
    simulation_years = collect(configured_simulation_years(config))
    start_year = first(simulation_years)
    end_year = last(simulation_years)
    source_years = management_source_years(config, simulation_years)
    nitrogen_deposition_year = nitrogen_deposition_source_year(config)
    landuse = read_management(
        catalog, :landuse, grid, cft_id;
        simulation_years = source_years, T = Float32, irrigated,
    )
    crop_mask = build_crop_mask(grid, landuse.values)
    phu_probe = read_management(
        catalog, :phu, grid, cft_id;
        selection = crop_mask.selection,
        simulation_years = source_years,
        active = falses(size(crop_mask.active)),
        T = Float32, irrigated,
    )
    valid_phu = map(axes(phu_probe.values, 2)) do cell
        active_years = view(crop_mask.active, :, cell)
        values = view(phu_probe.values, :, cell)
        all(!active || (isfinite(value) && value != 0) for (active, value) in zip(active_years, values))
    end
    any(valid_phu) || error("no landfrac-selected cells have valid PHU")
    excluded = .!valid_phu
    landfrac = vec(maximum(crop_mask.fraction; dims = 1))
    selection = select_cells(
        grid, crop_mask.selection.compact_indices[valid_phu],
    )
    cell_limit = Int(get(run, "cell_limit", 0))
    cell_limit >= 0 || error("cell_limit must be non-negative")
    if cell_limit > 0
        selection = select_cells(
            grid, selection.compact_indices[1:min(cell_limit, length(selection.cell_ids))],
        )
    end
    eligible_cell_count = length(selection.cell_ids)
    if partition_count != 1 || partition_rank != 0
        selection = partition_cell_selection(selection, partition_rank, partition_count)
    end
    excluded_landfrac = sum(landfrac[excluded])
    write_report(joinpath(output_directory, "selection_qc.toml"), Dict(
        "landfrac_positive_cells" => length(crop_mask.selection.cell_ids),
        "valid_phu_cells" => count(valid_phu),
        "excluded_phu_cells" => count(excluded),
        "excluded_landfrac_sum" => excluded_landfrac,
        "excluded_landfrac_fraction" => excluded_landfrac / sum(landfrac),
        "eligible_cells" => eligible_cell_count,
        "run_cells" => length(selection.cell_ids),
        "partition_rank" => Int(partition_rank),
        "partition_count" => Int(partition_count),
    ))
    hwsd = load_hwsd_targets(paths, selection)
    pool_allocation = haskey(paths, "pool_allocation") ?
        read_soil_pool_allocation(abspath(paths["pool_allocation"])) : nothing
    !isnothing(pool_allocation) && pool_allocation.cft_id != cft_id && error(
        "pool allocation CFT $(pool_allocation.cft_id) does not match requested CFT $cft_id",
    )
    !isnothing(pool_allocation) && pool_allocation.irrigated != irrigated && error(
        "pool allocation irrigation metadata does not match the requested water system",
    )
    selected_landuse = read_management(
        catalog, :landuse, grid, cft_id;
        selection, simulation_years = source_years, T = Float32, irrigated,
    )
    active = selected_landuse.values .> 0
    management = read_management_schedule(
        catalog, grid, selection, active, source_years, simulation_years, config, cft_id;
        irrigated,
    )
    management_blocks = [annual_management(management, index) for index in eachindex(simulation_years)]
    initial_data = model_inputs(
        grid, selection, hwsd, catalog, management; pool_allocation,
    )

    warmup_start_year = Int(get(run, "warmup_climate_start_year", start_year))
    warmup_end_year = Int(get(run, "warmup_climate_end_year", end_year))
    warmup_start_year <= warmup_end_year || error(
        "warmup_climate_start_year must not exceed warmup_climate_end_year",
    )
    warmup_years = collect(warmup_start_year:warmup_end_year)
    warmup_source_years = management_source_years(config, warmup_years)
    warmup_landuse = read_management(
        catalog, :landuse, grid, cft_id;
        selection, simulation_years = warmup_source_years, T = Float32, irrigated,
    )
    warmup_management = read_management_schedule(
        catalog, grid, selection, warmup_landuse.values .> 0,
        warmup_source_years, warmup_years, config, cft_id; irrigated,
    )
    warmup_management_blocks = [
        annual_management(warmup_management, index) for index in eachindex(warmup_years)
    ]
    warmup_reader = climate_blocks(
        catalog, grid;
        selection, start_year = warmup_start_year, end_year = warmup_end_year,
        block_days = 365, T = Float32, nitrogen_deposition_year,
    )
    climate_days(warmup_reader) == 365 * length(warmup_years) || error(
        "warm-up climate must contain complete 365-day years",
    )
    warmup_cache_climate = Bool(get(run, "warmup_cache_climate", false))
    warmup_forcings = climate_forcings(warmup_reader)
    if warmup_cache_climate
        println("caching $(length(warmup_forcings)) warm-up climate blocks in host memory")
        warmup_forcings = collect(warmup_forcings)
    end

    expected_days = production_enabled ? 365 * length(simulation_years) : 365
    reader = nothing
    forcings = nothing
    if production_enabled
        reader = climate_blocks(
            catalog, grid;
            selection, start_year, end_year, block_days = 365, T = Float32,
            nitrogen_deposition_year,
        )
        climate_days(reader) == expected_days || error("expected exactly $expected_days forcing days")
        first_block = read_climate_block(reader, 1)
        last_block = read_climate_block(reader, length(reader))
        start_date = first(first_block.time)
        end_date = last(last_block.time)
        (year(start_date), month(start_date), day(start_date)) == (start_year, 1, 1) ||
            error("forcing must start on $start_year-01-01, got $start_date")
        (year(end_date), month(end_date), day(end_date)) == (end_year, 12, 31) ||
            error("forcing must end on $end_year-12-31, got $end_date")
        forcings = climate_forcings(reader)
        length(forcings) == length(simulation_years) || error(
            "365-day blocks must produce one forcing block per simulation year",
        )
    end

    preflight_stream = OutputStream(
        production_output_variables();
        frequency = :annual, cell_ids = selection.cell_ids,
    )
    memory_safety_factor = Float64(get(run, "memory_safety_factor", 1.2))
    preflight_memory = estimate_memory(
        length(selection.cell_ids), expected_days;
        T = Float32,
        diagnostics = false,
        block_days = 365,
        backend = backend.name === :cpu ? :cpu : :accelerator,
        safety_factor = memory_safety_factor,
        warmup_years = Int(get(run, "warmup_maximum_years", get(run, "warmup_years", 10))),
        cached_forcing_blocks = warmup_cache_climate ? length(warmup_years) : 0,
        output_stream = preflight_stream,
    )
    write_report(joinpath(output_directory, "memory_preflight.toml"), Dict(
        string(name) => value for (name, value) in pairs(preflight_memory)
    ))
    if backend.name === :cuda
        available_device_bytes = CUDA.free_memory()
        preflight_memory.recommended_device_peak_bytes <= available_device_bytes || error(
            "estimated CUDA peak $(preflight_memory.recommended_device_peak_gib) GiB " *
            "exceeds available device memory $(available_device_bytes / 2.0^30) GiB",
        )
    end

    minimum_warmup_years = Int(get(run, "warmup_minimum_years", get(run, "warmup_years", 10)))
    maximum_warmup_years = Int(get(run, "warmup_maximum_years", minimum_warmup_years))
    warmup_options = (
        years = minimum_warmup_years,
        maximum_years = maximum_warmup_years,
        target_constrained = Bool(get(run, "warmup_target_constrained", true)),
        consecutive_years = Int(get(run, "warmup_consecutive_years", 3)),
        relative_tolerance = Float64(get(run, "warmup_relative_tolerance", 0.01)),
        pool_fraction_tolerance = Float64(get(
            run, "warmup_pool_fraction_tolerance", 0.01,
        )),
        required_converged_fraction = Float64(get(
            run, "warmup_required_converged_fraction", 1.0,
        )),
    )
    warmup_simulation = create_simulation(
        initial_data, selection, config, expected_days, backend.device, cft_id;
        irrigated, diagnostics = false,
    )
    warmup = agricultural_warmup!(
        warmup_simulation, warmup_forcings;
        warmup_options..., management_blocks = warmup_management_blocks,
        convergence_reducer = warmup_convergence_reducer,
    )
    warmup_drift = agricultural_warmup_drift(warmup)
    write_report(
        joinpath(output_directory, "warmup_cn_drift.toml"),
        warmup_drift,
    )
    write_warmup_cell_audit(
        joinpath(output_directory, "warmup_cell_audit.nc"), grid, selection, warmup,
    )
    write_soil_pool_allocation(
        joinpath(output_directory, "warmup_soil_pool_allocation.nc"),
        SoilPoolAllocation(
            selection,
            warmup.calibrated_pool_allocation.fast_carbon_fraction,
            warmup.calibrated_pool_allocation.fast_nitrogen_fraction,
            warmup.calibrated_pool_allocation.c_shift_fast,
            warmup.calibrated_pool_allocation.c_shift_slow;
            cft_id,
            irrigated,
            provenance = (
                source = "agricultural_warmup",
                warmup_years = warmup.years,
                target_constrained = warmup.target_constrained,
            ),
        ),
    )
    println(
        "warm-up completed: years=$(warmup.years), converged=$(warmup.converged), " *
        "converged_cell_fraction=$(warmup.converged_cell_fraction)",
    )
    require_warmup_convergence = Bool(get(run, "require_warmup_convergence", true))
    require_warmup_convergence && !warmup.converged && error(
        "warm-up did not reach the configured convergence requirement; " *
        "review warmup_cn_drift.toml or set require_warmup_convergence=false " *
        "for an explicitly non-production diagnostic run",
    )
    warmup_contract = (
        years = warmup.years,
        target_constrained = warmup.target_constrained,
        consecutive_years = warmup.consecutive_years,
        relative_tolerance = warmup.relative_tolerance,
        pool_fraction_tolerance = warmup.pool_fraction_tolerance,
        required_converged_fraction = warmup.required_converged_fraction,
        converged = warmup.converged,
    )
    if !production_enabled
        return (
            cells = length(selection.cell_ids),
            backend = backend.name,
            warmup_years = warmup_contract.years,
            warmup_converged = warmup_contract.converged,
            output_directory,
        )
    end
    warmup_checkpoint = joinpath(output_directory, "warmup_checkpoint.jld2")
    save_checkpoint(warmup_checkpoint, warmup_simulation)
    warmup = nothing
    warmup_drift = nothing
    GC.gc(true)
    backend.name === :cuda && CUDA.reclaim()

    production = create_simulation(
        initial_data, selection, config, expected_days, backend.device, cft_id;
        irrigated, diagnostics = false,
    )
    restore_checkpoint!(production, warmup_checkpoint)
    state_equal(production.state.prognostic, warmup_simulation.state.prognostic) ||
        error("warm-up checkpoint did not restore the prognostic state exactly")
    warmup_simulation = nothing
    GC.gc(true)
    backend.name === :cuda && CUDA.reclaim()

    annual_chunks = OutputChunk[]
    compact_writer = NetCDFBlockWriter(joinpath(output_directory, "compact"); prefix = "wheat")
    writer = chunk -> begin
        chunk.frequency === :annual && push!(annual_chunks, chunk)
        compact_writer(chunk)
    end
    stream = OutputStream(
        production_output_variables();
        frequency = :annual, writer, cell_ids = selection.cell_ids,
    )
    memory = estimate_memory(
        production, reader; prefetch = false, output_stream = stream,
    )
    write_report(joinpath(output_directory, "memory_estimate.toml"), Dict(
        string(name) => value for (name, value) in pairs(memory)
    ))

    run_output_block!(production, forcings[1], management_blocks[1], stream)
    production.simulated_days == 365 || error("$start_year run did not stop at day 365")
    first_year_checkpoint = joinpath(output_directory, "end_$(start_year)_checkpoint.jld2")
    save_checkpoint(first_year_checkpoint, production)

    continued = create_simulation(
        initial_data, selection, config, expected_days, backend.device, cft_id;
        irrigated, diagnostics = false,
    )
    restore_checkpoint!(continued, first_year_checkpoint)
    state_equal(continued.state.prognostic, production.state.prognostic) ||
        error("$start_year checkpoint did not restore the prognostic state exactly")
    for index in 2:length(forcings)
        run_output_block!(continued, forcings[index], management_blocks[index], stream)
    end
    finish_output_stream!(stream, continued.simulated_days)
    continued.simulated_days == expected_days ||
        error("production run did not reach the end of $end_year")
    soil_water = Array(continued.state.prognostic.soil.water.storage)
    all(isfinite, soil_water) || error("non-finite soil water")
    all(>=(0), soil_water) || error("negative soil water")
    for (name, pool) in pairs(continued.state.prognostic.soil.carbon)
        host_pool = Array(pool)
        invalid = findfirst(value -> !isfinite(value) || value < 0, host_pool)
        if invalid !== nothing
            compact_cell = Tuple(invalid)[end]
            error(
                "invalid soil carbon pool: pool=$(name), index=$(invalid), " *
                "cell_id=$(selection.cell_ids[compact_cell]), value=$(host_pool[invalid])",
            )
        end
    end
    for (name, pool) in pairs(continued.state.prognostic.soil.nitrogen)
        host_pool = Array(pool)
        invalid = findfirst(value -> !isfinite(value) || value < 0, host_pool)
        if invalid !== nothing
            compact_cell = Tuple(invalid)[end]
            error(
                "invalid soil nitrogen pool: pool=$(name), index=$(invalid), " *
                "cell_id=$(selection.cell_ids[compact_cell]), value=$(host_pool[invalid])",
            )
        end
    end
    write_reconstructed_output(
        joinpath(output_directory, "global_wheat_$(start_year)_$(end_year).nc"),
        grid, selection, annual_chunks, simulation_years,
    )

    diagnostic_count = min(Int(run["diagnostic_cells"]), length(selection.cell_ids))
    diagnostic_selection = select_cells(grid, selection.compact_indices[1:diagnostic_count])
    diagnostic_landuse = read_management(
        catalog, :landuse, grid, cft_id;
        selection = diagnostic_selection,
        simulation_years = source_years,
        T = Float32, irrigated,
    )
    diagnostic_management = read_management_schedule(
        catalog, grid, diagnostic_selection, diagnostic_landuse.values .> 0,
        source_years, simulation_years, config, cft_id; irrigated,
    )
    diagnostic_management_blocks = [
        annual_management(diagnostic_management, index) for index in eachindex(simulation_years)
    ]
    diagnostic_initial = model_inputs(
        grid, diagnostic_selection, hwsd, catalog, diagnostic_management;
        pool_allocation,
    )
    diagnostic_reader = climate_blocks(
        catalog, grid;
        selection = diagnostic_selection, start_year, end_year,
        block_days = 365, T = Float32, nitrogen_deposition_year,
    )
    diagnostic = create_simulation(
        diagnostic_initial, diagnostic_selection, config, expected_days, backend.device, cft_id;
        irrigated, diagnostics = !irrigated,
    )
    diagnostic_warmup_reader = climate_blocks(
        catalog, grid;
        selection = diagnostic_selection,
        start_year = warmup_start_year, end_year = warmup_end_year,
        block_days = 365, T = Float32, nitrogen_deposition_year,
    )
    diagnostic_warmup_landuse = read_management(
        catalog, :landuse, grid, cft_id;
        selection = diagnostic_selection,
        simulation_years = warmup_source_years,
        T = Float32, irrigated,
    )
    diagnostic_warmup_management = read_management_schedule(
        catalog, grid, diagnostic_selection,
        diagnostic_warmup_landuse.values .> 0,
        warmup_source_years, warmup_years, config, cft_id; irrigated,
    )
    diagnostic_warmup_management_blocks = [
        annual_management(diagnostic_warmup_management, index)
        for index in eachindex(warmup_years)
    ]
    agricultural_warmup!(
        diagnostic, climate_forcings(diagnostic_warmup_reader);
        years = warmup_contract.years,
        maximum_years = warmup_contract.years,
        target_constrained = warmup_contract.target_constrained,
        consecutive_years = warmup_contract.consecutive_years,
        relative_tolerance = warmup_contract.relative_tolerance,
        pool_fraction_tolerance = warmup_contract.pool_fraction_tolerance,
        required_converged_fraction = warmup_contract.required_converged_fraction,
        management_blocks = diagnostic_warmup_management_blocks,
    )
    run_simulation!(
        diagnostic, climate_forcings(diagnostic_reader);
        spinup = false, management_blocks = diagnostic_management_blocks,
    )
    balance = if irrigated
        Dict{String, Any}(
            "water_balance" => "not evaluated for irrigated simulations",
            "diagnostic_summary" => "daily balance ledgers are disabled for irrigation",
        )
    else
        report = balance_report(diagnostic, config["validation"])
        report["simulation_summary"] = simulation_summary(diagnostic)
        report
    end
    write_report(joinpath(output_directory, "sampled_balance_summary.toml"), balance)
    return (;
        cells = length(selection.cell_ids),
        backend = backend.name,
        warmup_years = warmup_contract.years,
        warmup_converged = warmup_contract.converged,
        memory,
        output_directory,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: run_global_wheat_cpu.jl CONFIG_TOML")
    println(run_global_wheat(abspath(ARGS[1]); backend_override = :cpu))
end
