using MPI

include(joinpath(@__DIR__, "run_global_cfts_cpu.jl"))

function mpi_warmup_convergence_reducer(comm, converged_cells, total_cells)
    global_converged = MPI.Allreduce(Int64(converged_cells), +, comm)
    global_total = MPI.Allreduce(Int64(total_cells), +, comm)
    return global_converged, global_total
end

function mpi_resume_consensus(comm, locally_ready)
    ready_ranks = MPI.Allreduce(locally_ready ? Int32(1) : Int32(0), +, comm)
    return ready_ranks == MPI.Comm_size(comm)
end

function merge_partition_allocations(paths, output_path)
    allocations = read_soil_pool_allocation.(paths)
    cft_id = only(unique(getproperty.(allocations, :cft_id)))
    irrigated = only(unique(getproperty.(allocations, :irrigated)))
    cell_ids = reduce(vcat, getproperty.(getproperty.(allocations, :selection), :cell_ids))
    allunique(cell_ids) || error("MPI allocation partitions contain overlapping cells")
    order = sortperm(cell_ids)
    values(name) = reduce(hcat, getproperty.(allocations, name))[:, order]
    selection = NeuralCropData.CellSelection(collect(eachindex(cell_ids)), cell_ids[order])
    allocation = SoilPoolAllocation(
        selection,
        values(:fast_carbon_fraction),
        values(:fast_nitrogen_fraction),
        values(:c_shift_fast),
        values(:c_shift_slow);
        cft_id,
        irrigated,
        provenance = (
            source = "mpi_partition_merge",
            partition_count = length(paths),
        ),
    )
    mkpath(dirname(output_path))
    return write_soil_pool_allocation(output_path, allocation)
end

function merge_partition_outputs(output_paths, allocation_paths, output_path)
    length(output_paths) == length(allocation_paths) || throw(DimensionMismatch(
        "MPI output and allocation partition counts differ",
    ))
    allocations = read_soil_pool_allocation.(allocation_paths)
    reference = NCDataset(first(output_paths), "r")
    try
        longitude = Float32.(reference["longitude"][:])
        latitude = Float32.(reference["latitude"][:])
        time = Int32.(reference["time"][:])
        cellid = Int32.(reference["cellid"][:, :])
        coordinate_names = Set(("longitude", "latitude", "time", "cellid"))
        variable_names = sort!(filter(
            name -> !(name in coordinate_names),
            String.(collect(keys(reference))),
        ))
        expected_names = union(coordinate_names, Set(variable_names))
        merged = Dict(name => fill(
            Float32(NaN), length(longitude), length(latitude), length(time),
        ) for name in variable_names)
        positions = Dict{Int32, CartesianIndex{2}}()
        for index in CartesianIndices(cellid)
            positions[cellid[index]] = index
        end
        assigned = Set{Int32}()

        for (path, allocation) in zip(output_paths, allocations)
            NCDataset(path, "r") do dataset
                Set(String.(collect(keys(dataset)))) == expected_names || error(
                    "output variables differ across MPI partitions",
                )
                dataset["longitude"][:] == longitude || error(
                    "longitude coordinates differ across MPI outputs",
                )
                dataset["latitude"][:] == latitude || error(
                    "latitude coordinates differ across MPI outputs",
                )
                Int32.(dataset["time"][:]) == time || error(
                    "time coordinates differ across MPI outputs",
                )
                Int32.(dataset["cellid"][:, :]) == cellid || error(
                    "cell IDs differ across MPI outputs",
                )
                indices = map(allocation.selection.cell_ids) do cell_id
                    cell_id in assigned && error("MPI output partitions overlap at cell $cell_id")
                    haskey(positions, cell_id) || error(
                        "MPI output cell $cell_id is absent from the canonical grid",
                    )
                    push!(assigned, cell_id)
                    return positions[cell_id]
                end
                for name in variable_names
                    size(dataset[name]) == size(merged[name]) || error(
                        "variable $name dimensions differ across MPI outputs",
                    )
                    partition_values = dataset[name][:, :, :]
                    for index in indices
                        merged[name][index[1], index[2], :] .=
                            partition_values[index[1], index[2], :]
                    end
                end
            end
        end

        mkpath(dirname(output_path))
        isfile(output_path) && rm(output_path; force = true)
        NCDataset(output_path, "c") do dataset
            defDim(dataset, "longitude", length(longitude))
            defDim(dataset, "latitude", length(latitude))
            defDim(dataset, "time", length(time))
            defVar(dataset, "longitude", Float32, ("longitude",))[:] = longitude
            defVar(dataset, "latitude", Float32, ("latitude",))[:] = latitude
            defVar(dataset, "time", Int32, ("time",))[:] = time
            defVar(dataset, "cellid", Int32, ("longitude", "latitude"))[:, :] = cellid
            for name in variable_names
                defVar(
                    dataset, name, Float32, ("longitude", "latitude", "time"),
                )[:, :, :] = merged[name]
            end
            dataset.attrib["source"] = "mpi_partition_merge"
            dataset.attrib["partition_count"] = length(output_paths)
        end
    finally
        close(reference)
    end
    return output_path
end

function merge_mpi_rank_products(
    rank_manifests,
    output_root;
    manifest_path = joinpath(output_root, "cft_batch_manifest.toml"),
)
    manifests = TOML.parsefile.(rank_manifests)
    batch_maps = [Dict(String(batch["name"]) => batch for batch in manifest["batches"])
                  for manifest in manifests]
    names = sort!(collect(keys(first(batch_maps))))
    all(Set(keys(batch_map)) == Set(names) for batch_map in batch_maps) || error(
        "MPI ranks completed different CFT batches",
    )
    products = Dict{String, Any}[]
    for name in names
        batches = [batch_map[name] for batch_map in batch_maps]
        all(batch["cft_id"] == first(batches)["cft_id"] for batch in batches) || error(
            "CFT IDs differ across MPI rank manifests for $name",
        )
        all(batch["irrigated"] == first(batches)["irrigated"] for batch in batches) || error(
            "water systems differ across MPI rank manifests for $name",
        )
        allocation_paths = String.(getindex.(batches, "pool_allocation"))
        all(isfile, allocation_paths) || error("MPI allocation partition missing for $name")
        # Publish merged MPI products through the same batches/<patch> contract
        # as the single-device runner. Rank-local products remain temporary.
        merged_directory = joinpath(output_root, "batches", name)
        merged_allocation = merge_partition_allocations(
            allocation_paths,
            joinpath(
                merged_directory, "calibration", "warmup_soil_pool_allocation.nc",
            ),
        )
        output_paths = String.(getindex.(batches, "production_output"))
        output_exists = isfile.(output_paths)
        any(output_exists) && !all(output_exists) && error(
            "only some MPI production partitions exist for $name",
        )
        merged_output = if all(output_exists)
            merge_partition_outputs(
                output_paths,
                allocation_paths,
                joinpath(
                    merged_directory, "production", basename(first(output_paths)),
                ),
            )
        else
            ""
        end
        push!(products, Dict(
            "name" => name,
            "cft_id" => first(batches)["cft_id"],
            "irrigated" => first(batches)["irrigated"],
            "pool_allocation" => abspath(merged_allocation),
            "production_output" => isempty(merged_output) ? "" : abspath(merged_output),
        ))
    end
    write_report(manifest_path, Dict(
        "schema_version" => "1",
        "repository_commit" => _repository_commit(),
        "partition_count" => length(rank_manifests),
        "batches" => products,
        "created_at" => string(now()),
    ))
    return manifest_path
end

function run_global_cfts_mpi(
    config_path;
    cft_id::Union{Nothing, Integer} = nothing,
    irrigated::Union{Nothing, Bool} = nothing,
)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    count = MPI.Comm_size(comm)
    config = TOML.parsefile(config_path)
    configured_output = abspath(config["paths"]["output_directory"])
    mpi_config = get(config, "mpi", Dict{String, Any}())
    output_base = abspath(get(mpi_config, "output_directory", configured_output))
    cleanup_rank_outputs = Bool(get(mpi_config, "cleanup_rank_outputs", false))
    scoped_output = isnothing(cft_id) ? output_base :
        joinpath(output_base, "batches", batch_name(cft_id, irrigated))
    output_root = joinpath(scoped_output, "mpi_$(lpad(count, 4, '0'))_ranks")
    mpi_manifest_path = joinpath(output_root, "mpi_manifest.toml")
    # Independent CFT jobs may run concurrently. Keep the resume manifest
    # patch-local so jobs never overwrite one another.
    merged_manifest = joinpath(scoped_output, "cft_batch_manifest.toml")
    if cleanup_rank_outputs
        locally_complete = false
        if isfile(mpi_manifest_path)
            manifest = TOML.parsefile(mpi_manifest_path)
            merged_manifest = String(get(manifest, "merged_manifest", ""))
            locally_complete = get(manifest, "config_fingerprint", "") ==
                               _file_fingerprint(config_path) &&
                               get(manifest, "rank_count", 0) == count &&
                               isfile(merged_manifest)
        end
        if mpi_resume_consensus(comm, locally_complete)
            rank == 0 && println("MPI batch already merged: $merged_manifest")
            return merged_manifest
        end
    end
    rank_directory = joinpath(output_root, "rank_$(lpad(rank, 4, '0'))")
    convergence_reducer = (converged_cells, total_cells) ->
        mpi_warmup_convergence_reducer(comm, converged_cells, total_cells)
    resume_consensus = locally_ready -> mpi_resume_consensus(comm, locally_ready)

    println("MPI rank $rank/$count: output=$rank_directory")
    rank_manifest = run_global_cfts(
        config_path;
        backend_override = :cpu,
        cft_id,
        irrigated,
        output_directory_override = rank_directory,
        partition_rank = rank,
        partition_count = count,
        warmup_convergence_reducer = convergence_reducer,
        resume_consensus,
    )
    MPI.Barrier(comm)

    if rank == 0
        rank_directories = [joinpath(output_root, "rank_$(lpad(value, 4, '0'))")
                            for value in 0:(count - 1)]
        rank_manifests = [joinpath(path, "cft_batch_manifest.toml")
                          for path in rank_directories]
        all(isfile, rank_manifests) || error("one or more MPI rank manifests are missing")
        merge_mpi_rank_products(
            rank_manifests,
            output_base;
            manifest_path = merged_manifest,
        ) == merged_manifest || error(
            "unexpected merged MPI manifest path",
        )
        write_report(mpi_manifest_path, Dict(
            "schema_version" => "1",
            "repository_commit" => _repository_commit(),
            "config_path" => abspath(config_path),
            "config_fingerprint" => _file_fingerprint(config_path),
            "rank_count" => count,
            "rank_outputs_retained" => !cleanup_rank_outputs,
            "rank_directories" => rank_directories,
            "rank_manifests" => rank_manifests,
            "merged_manifest" => merged_manifest,
            "created_at" => string(now()),
        ))
        if cleanup_rank_outputs
            foreach(rank_directories) do directory
                rm(directory; recursive = true, force = true)
            end
        end
    end
    MPI.Barrier(comm)
    return merged_manifest
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (1, 3) || error(
        "usage: run_global_cfts_mpi.jl CONFIG_TOML [CFT_ID rainfed|irrigated]",
    )
    cft_id = length(ARGS) == 3 ? parse(Int, ARGS[2]) : nothing
    irrigated = length(ARGS) == 3 ?
        Symbol(lowercase(ARGS[3])) === :irrigated : nothing
    length(ARGS) == 3 && Symbol(lowercase(ARGS[3])) in (:rainfed, :irrigated) ||
        length(ARGS) == 1 || error("water system must be rainfed or irrigated")

    MPI.Init()
    try
        result = run_global_cfts_mpi(abspath(ARGS[1]); cft_id, irrigated)
        MPI.Comm_rank(MPI.COMM_WORLD) == 0 && println(result)
    catch exception
        rank = MPI.Comm_rank(MPI.COMM_WORLD)
        println(stderr, "MPI rank $rank failed: ", sprint(showerror, exception, catch_backtrace()))
        MPI.Abort(MPI.COMM_WORLD, 1)
    finally
        MPI.Finalize()
    end
end
