include(joinpath(@__DIR__, "..", "..", "scripts", "run_global_cfts_cpu.jl"))

@testset "multi-CFT patch selection, state isolation, and batch configuration" begin
    high_throughput_closure = percolation_energy_closure(
        Float32[10], Float32[2.0e6], Float32[0], Float32[0], Float32[0];
        absolute_tolerance = 5.0f0, relative_tolerance = 5.0f-6,
    )
    @test high_throughput_closure.closes
    @test high_throughput_closure.maximum_relative_residual ≈ 5.0f-6
    @test !percolation_energy_closure(
        Float32[20], Float32[2.0e6], Float32[0], Float32[0], Float32[0];
        absolute_tolerance = 5.0f0, relative_tolerance = 5.0f-6,
    ).closes

    @test requested_crop_systems(Dict{String, Any}()) == [(1, false)]
    @test requested_crop_systems(Dict("cfts" => Dict(
        "cft_ids" => "all", "water_systems" => ["rainfed", "irrigated"],
    ))) == [(cft_id, irrigated) for cft_id in 1:length(CFTS) for irrigated in (false, true)]
    selection = NeuralCropData.CellSelection(11:20, 101:110)
    partitions = [partition_cell_selection(selection, rank, 3) for rank in 0:2]
    @test length.(getproperty.(partitions, :cell_ids)) == [4, 3, 3]
    @test reduce(vcat, getproperty.(partitions, :cell_ids)) == selection.cell_ids
    @test all(isempty(intersect(partitions[i].cell_ids, partitions[j].cell_ids))
        for i in eachindex(partitions) for j in (i + 1):length(partitions))
    @test_throws ArgumentError partition_cell_selection(selection, 3, 3)
    @test_throws ArgumentError partition_cell_selection(selection, 0, 11)
    patch_domain = combine_patch_domains([
        PatchDomain([1], [3], [42], [1], [false], Float32[0.2]),
        PatchDomain([1], [3], [42], [2], [true], Float32[0.3]),
    ])
    @test patch_domain.cell_ids == Int32[42, 42]
    @test patch_domain.cft_ids == Int32[1, 2]
    @test patch_domain.irrigated == BitVector([false, true])

    single_catalog = catalog_from_config(Dict("paths" => Dict("input_directory" => "/tmp/input")))
    @test single_catalog.cfts.ids == Int32.(1:12)
    @test endswith(dataset(single_catalog, :landuse).path, "landuse_24cfts_2015.nc")
    @test cft_index(single_catalog.cfts, 1) == 1

    historical_catalog = catalog_from_config(Dict(
        "paths" => Dict("input_directory" => "/tmp/input"),
        "management" => Dict(
            "landuse_file" => "landuse_24cfts_1500-2017.nc",
            "fertilizer_file" => "fertilizer_24cfts_1860-2017.nc",
            "manure_file" => "manure_24cfts_1860-2017.nc",
            "residue_fraction_file" => "residue_12cfts_1850-2015.nc",
            "sowing_date_file" => "sdate_24cfts_static.nc",
            "phu_file" => "phu_24cfts_1901-2019.nc",
        ),
    ))
    @test basename(dataset(historical_catalog, :landuse).path) ==
        "landuse_24cfts_1500-2017.nc"
    @test basename(dataset(historical_catalog, :fertilizer).path) ==
        "fertilizer_24cfts_1860-2017.nc"
    @test basename(dataset(historical_catalog, :manure).path) ==
        "manure_24cfts_1860-2017.nc"
    @test basename(dataset(historical_catalog, :residue_fraction).path) ==
        "residue_12cfts_1850-2015.nc"
    @test basename(dataset(historical_catalog, :sowing_date).path) ==
        "sdate_24cfts_static.nc"
    @test basename(dataset(historical_catalog, :phu).path) ==
        "phu_24cfts_1901-2019.nc"

    days = 30
    rainfed = initialize_simulation(
        cft1, lifecycle_initial_data(Float32);
        days, diagnostics = false, irrigation = false, fertilizer = :yes,
    )
    irrigated = initialize_simulation(
        cft1, lifecycle_initial_data(Float32);
        days, diagnostics = false, irrigation = true, fertilizer = :yes,
    )
    @test rainfed.state.prognostic.soil.water.storage !==
        irrigated.state.prognostic.soil.water.storage
    initial_irrigated_water = copy(irrigated.state.prognostic.soil.water.storage)
    run_simulation!(rainfed, lifecycle_climate(Float32, days); spinup = false)
    @test irrigated.state.prognostic.soil.water.storage == initial_irrigated_water
    run_simulation!(irrigated, lifecycle_climate(Float32, days); spinup = false)

    mktempdir() do directory
        batch_config = Dict{String, Any}(
            "paths" => Dict{String, Any}("output_directory" => directory),
            "run" => Dict{String, Any}(
                "warmup_minimum_years" => 1,
                "warmup_maximum_years" => 2,
                "require_warmup_convergence" => true,
            ),
        )
        calibration_config = write_batch_config(
            joinpath(directory, "calibration.toml"), batch_config;
            output_directory = joinpath(directory, "calibration"),
            production = false,
            target_constrained = true,
        )
        calibration_run = TOML.parsefile(calibration_config)["run"]
        @test calibration_run["warmup_minimum_years"] == 600
        @test calibration_run["warmup_maximum_years"] == 600
        @test !calibration_run["require_warmup_convergence"]

        free_config = write_batch_config(
            joinpath(directory, "production.toml"), batch_config;
            output_directory = joinpath(directory, "production"),
            pool_allocation = joinpath(directory, "allocation.nc"),
            production = true,
            target_constrained = false,
            warmup_options = Dict{String, Any}(
                "warmup_minimum_years" => 600,
                "warmup_maximum_years" => 1500,
                "require_warmup_convergence" => true,
            ),
        )
        free_run = TOML.parsefile(free_config)["run"]
        @test free_run["warmup_minimum_years"] == 600
        @test free_run["warmup_maximum_years"] == 1500
        @test free_run["require_warmup_convergence"]

        allocation_path = joinpath(directory, "allocation.nc")
        production_output = joinpath(directory, "production.nc")
        write(allocation_path, "allocation")
        write(production_output, "production")
        status_path = joinpath(directory, "batch_status.toml")
        fingerprint = _file_fingerprint(free_config)
        _write_batch_status(
            status_path;
            status = "calibration_completed",
            name = "cft_01_rainfed",
            cft_id = 1,
            irrigated = false,
            calibration_config,
            production_config = free_config,
            production_output,
            calibration_output_directory = joinpath(directory, "calibration"),
            production_output_directory = joinpath(directory, "production"),
            allocation_path,
            production_config_fingerprint = fingerprint,
            allocation_fingerprint = _file_fingerprint(allocation_path),
            calibration_completed = true,
        )
        @test _resumable_calibration(status_path, fingerprint) !== nothing
        @test _resumable_batch(status_path, fingerprint) === nothing

        _write_batch_status(
            status_path;
            status = "failed",
            name = "cft_01_rainfed",
            cft_id = 1,
            irrigated = false,
            calibration_config,
            production_config = free_config,
            production_output,
            calibration_output_directory = joinpath(directory, "calibration"),
            production_output_directory = joinpath(directory, "production"),
            allocation_path,
            production_config_fingerprint = fingerprint,
            allocation_fingerprint = _file_fingerprint(allocation_path),
            calibration_completed = true,
        )
        @test _resumable_calibration(status_path, fingerprint) !== nothing

        _write_batch_status(
            status_path;
            status = "completed",
            name = "cft_01_rainfed",
            cft_id = 1,
            irrigated = false,
            calibration_config,
            production_config = free_config,
            production_output,
            calibration_output_directory = joinpath(directory, "calibration"),
            production_output_directory = joinpath(directory, "production"),
            allocation_path,
            production_config_fingerprint = fingerprint,
        )
        @test _resumable_batch(status_path, fingerprint) !== nothing
        @test _resumable_batch(status_path, "different-config") === nothing
    end

end
