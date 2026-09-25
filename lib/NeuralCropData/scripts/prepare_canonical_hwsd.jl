push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))

using NeuralCropData
using NCDatasets
using TOML

include(joinpath(@__DIR__, "extract_hwsd_cell.jl"))

function query_all_layer_rows(data_directory)
    csv_path = joinpath(data_directory, "HWSD2_LAYERS.csv")
    if isfile(csv_path)
        lines = filter(!isempty, strip.(readlines(csv_path)))
        isempty(lines) && error("HWSD layer CSV contains no rows: $csv_path")
        return [split(line, ','; keepempty = true) for line in lines]
    end

    database_path = joinpath(data_directory, "HWSD2.mdb")
    mdb_sql = Sys.which("mdb-sql")
    isnothing(mdb_sql) && error(
        "mdb-sql is unavailable and $csv_path was not found. " *
        "Export HWSD2_LAYERS.csv on a machine with MDB Tools and copy it " *
        "into the HWSD directory.",
    )
    sql = """
        SELECT HWSD2_SMU_ID,SEQUENCE,SHARE,ROOT_DEPTH,WRB4,LAYER,COARSE,BULK,ORG_CARBON,TOTAL_N
        FROM HWSD2_LAYERS;
        """
    command = pipeline(
        `$mdb_sql -P -H -F -d , $database_path`,
        stdin = IOBuffer(sql),
    )
    lines = filter(!isempty, strip.(readlines(command)))
    isempty(lines) && error("HWSD layer query returned no rows")
    return [split(line, ','; keepempty = true) for line in lines]
end

function indexed_mapping_units(targets_by_unit)
    unit_count = Int(SOURCE_NODATA) + 1
    carbon = fill(NaN, 5, unit_count)
    nitrogen = fill(NaN, 5, unit_count)
    uncertain = falses(5, unit_count)
    soil = falses(unit_count)
    for (mapping_unit, target) in targets_by_unit
        index = mapping_unit + 1
        checkbounds(Bool, soil, index) || continue
        soil[index] = target.is_soil
        carbon[:, index] = target.carbon
        nitrogen[:, index] = target.nitrogen
        uncertain[:, index] = target.uncertain
    end
    return (; carbon, nitrogen, uncertain, soil)
end

function accumulate_canonical_grid!(accumulator, raster_path, grid, lookup)
    cells_by_latitude = [Int[] for _ in eachindex(grid.latitude)]
    for compact_index in eachindex(grid.cell_ids)
        push!(cells_by_latitude[Int(grid.latitude_indices[compact_index])], compact_index)
    end

    row_values = Vector{UInt16}(undef, SOURCE_COLUMNS)
    longitude_starts = floor.(
        Int,
        (grid.longitude[Int.(grid.longitude_indices)] .- 0.25 .+ 180) ./ SOURCE_RESOLUTION,
    ) .+ 1
    open(raster_path, "r") do io
        for latitude_index in eachindex(grid.latitude)
            isempty(cells_by_latitude[latitude_index]) && continue
            south = grid.latitude[latitude_index] - 0.25
            for source_latitude_index in 1:60
                source_latitude = south + (source_latitude_index - 0.5) * SOURCE_RESOLUTION
                source_row = floor(Int, (90 - source_latitude) / SOURCE_RESOLUTION) + 1
                seek(io, (source_row - 1) * SOURCE_ROW_BYTES)
                read!(io, row_values)
                pixel_area = EARTH_RADIUS_M^2 * deg2rad(SOURCE_RESOLUTION) * (
                    sin(deg2rad(source_latitude + SOURCE_RESOLUTION / 2)) -
                    sin(deg2rad(source_latitude - SOURCE_RESOLUTION / 2))
                )

                for cell in cells_by_latitude[latitude_index]
                    first_column = longitude_starts[cell]
                    for column in first_column:(first_column + 59)
                        mapping_unit = Int(ltoh(row_values[column]))
                        unit_index = mapping_unit + 1
                        lookup.soil[unit_index] || continue
                        accumulator.target_area[cell] += pixel_area
                        for layer in 1:5
                            carbon = lookup.carbon[layer, unit_index]
                            if isfinite(carbon)
                                accumulator.carbon_sum[layer, cell] += carbon * pixel_area
                                accumulator.carbon_area[layer, cell] += pixel_area
                            end
                            nitrogen = lookup.nitrogen[layer, unit_index]
                            if isfinite(nitrogen)
                                accumulator.nitrogen_sum[layer, cell] += nitrogen * pixel_area
                                accumulator.nitrogen_area[layer, cell] += pixel_area
                            end
                            lookup.uncertain[layer, unit_index] &&
                                (accumulator.uncertain_area[layer, cell] += pixel_area)
                        end
                    end
                end
            end
            latitude_index % 10 == 0 && println(
                "processed latitude band $latitude_index / $(length(grid.latitude))",
            )
        end
    end
    return accumulator
end

function nearest_complete_cell(grid, complete, cell, max_radius_cells)
    longitude_index = Int(grid.longitude_indices[cell])
    latitude_index = Int(grid.latitude_indices[cell])
    best = nothing
    for radius in 1:max_radius_cells
        for longitude_offset in -radius:radius, latitude_offset in -radius:radius
            max(abs(longitude_offset), abs(latitude_offset)) == radius || continue
            candidate_longitude = mod1(longitude_index + longitude_offset, length(grid.longitude))
            candidate_latitude = latitude_index + latitude_offset
            checkbounds(Bool, grid.latitude, candidate_latitude) || continue
            cell_id = grid.cellid[candidate_longitude, candidate_latitude]
            cell_id >= 0 || continue
            candidate = searchsortedfirst(grid.cell_ids, cell_id)
            checkbounds(Bool, complete, candidate) && complete[candidate] || continue
            distance = distance_km(
                grid.longitude[longitude_index], grid.latitude[latitude_index],
                grid.longitude[candidate_longitude], grid.latitude[candidate_latitude],
            )
            if isnothing(best) || distance < best.distance
                best = (; cell = candidate, distance)
            end
        end
        !isnothing(best) && return best
    end
    return nothing
end

function apply_fallback(targets, grid; max_radius_cells, production_active)
    carbon = copy(targets.soil_organic_carbon)
    nitrogen = copy(targets.total_nitrogen)
    coverage = copy(targets.coverage)
    uncertain = copy(targets.uncertain)
    complete = vec(all(isfinite.(carbon) .& isfinite.(nitrogen); dims = 1))
    length(production_active) == length(complete) ||
        throw(DimensionMismatch("production mask must match canonical cells"))
    fallback_used = falses(length(complete))
    donor_cell_id = fill(Int32(-1), length(complete))
    fill_distance_km = fill(Float32(NaN), length(complete))
    original_minimum_coverage = Float32.(vec(minimum(coverage; dims = 1)))

    for cell in eachindex(complete)
        complete[cell] && continue
        donor = nearest_complete_cell(grid, complete, cell, max_radius_cells)
        isnothing(donor) && continue
        carbon[:, cell] = carbon[:, donor.cell]
        nitrogen[:, cell] = nitrogen[:, donor.cell]
        coverage[:, cell] = coverage[:, donor.cell]
        uncertain[:, cell] .= true
        fallback_used[cell] = true
        donor_cell_id[cell] = grid.cell_ids[donor.cell]
        fill_distance_km[cell] = donor.distance
    end
    unresolved = .!complete .& .!fallback_used
    provenance = merge(targets.provenance, (
        fill_policy = "nearest complete canonical 0.5-degree HWSD cell",
    ))
    filled = SoilCNTargets(
        targets.selection, targets.layer_bounds, carbon, nitrogen, coverage,
        uncertain, provenance,
    )
    qc = (;
        complete_before_fallback = count(complete),
        fallback_cells = count(fallback_used),
        unresolved_cells = count(unresolved),
        production_active_cells = count(production_active),
        active_fallback_cells = count(fallback_used .& production_active),
        active_unresolved_cells = count(unresolved .& production_active),
        inactive_unresolved_cells = count(unresolved .& .!production_active),
        fallback_used,
        donor_cell_id,
        fill_distance_km,
        original_minimum_coverage,
    )
    return filled, qc
end

relative_error(actual, reference) = reference == 0 ? abs(actual) : abs(actual - reference) / abs(reference)

function write_qc(
    path, output_path, grid, accumulator, targets, fallback, production_mask_source,
)
    complete = .!fallback.fallback_used .& vec(all(isfinite.(
        targets.soil_organic_carbon,
    ) .& isfinite.(targets.total_nitrogen); dims = 1))
    source_carbon = vec(sum(accumulator.carbon_sum[:, complete]; dims = 2))
    source_nitrogen = vec(sum(accumulator.nitrogen_sum[:, complete]; dims = 2))
    output_carbon = [sum(
        targets.soil_organic_carbon[layer, complete] .*
        accumulator.carbon_area[layer, complete],
    ) for layer in 1:5]
    output_nitrogen = [sum(
        targets.total_nitrogen[layer, complete] .*
        accumulator.nitrogen_area[layer, complete],
    ) for layer in 1:5]
    finite_distances = filter(isfinite, fallback.fill_distance_km)
    report = Dict(
        "schema_version" => string(DATA_SCHEMA_VERSION),
        "source_version" => "HWSD v2.01",
        "output" => abspath(output_path),
        "production_mask_source" => production_mask_source,
        "canonical_cells" => length(grid.cell_ids),
        "complete_before_fallback" => fallback.complete_before_fallback,
        "fallback_cells" => fallback.fallback_cells,
        "unresolved_cells" => fallback.unresolved_cells,
        "production_active_cells" => fallback.production_active_cells,
        "active_fallback_cells" => fallback.active_fallback_cells,
        "active_unresolved_cells" => fallback.active_unresolved_cells,
        "inactive_unresolved_cells" => fallback.inactive_unresolved_cells,
        "minimum_coverage" => minimum(targets.coverage),
        "uncertain_layer_cells" => count(targets.uncertain),
        "maximum_fallback_distance_km" =>
            isempty(finite_distances) ? 0.0f0 : maximum(finite_distances),
        "carbon_source_g" => source_carbon,
        "nitrogen_source_g" => source_nitrogen,
        "carbon_output_g" => output_carbon,
        "nitrogen_output_g" => output_nitrogen,
        "carbon_relative_error" => relative_error.(output_carbon, source_carbon),
        "nitrogen_relative_error" => relative_error.(output_nitrogen, source_nitrogen),
    )
    open(path, "w") do io
        TOML.print(io, report; sorted = true)
    end
    return report
end

function prepare_canonical_hwsd(
    data_directory, grid_path, output_path, qc_path;
    max_fill_radius_cells = 8,
    landuse_path = nothing,
)
    grid = read_grid(grid_path; T = Float64)
    production_mask_source = isnothing(landuse_path) ?
        "all canonical cells" : abspath(landuse_path)
    production_active = if isnothing(landuse_path)
        trues(length(grid.cell_ids))
    else
        spec = DatasetSpec(landuse_path, "landfrac"; cft_ids = [1])
        registry = CFTRegistry([1], ["rainfed wheat"])
        landuse = read_management(spec, :landuse, grid, registry, 1; T = Float32)
        size(landuse.values, 1) == 1 || error(
            "production land-use file must contain exactly one selected year",
        )
        BitVector(vec(landuse.values .> 0))
    end
    mkpath(dirname(output_path))
    mkpath(dirname(qc_path))
    rows = query_all_layer_rows(data_directory)
    lookup = indexed_mapping_units(mapping_unit_targets(rows))
    accumulator = init_soil_cn_aggregator(Float64, 5, length(grid.cell_ids))
    accumulate_canonical_grid!(
        accumulator, joinpath(data_directory, "HWSD2.bil"), grid, lookup,
    )
    provenance = (
        schema_version = DATA_SCHEMA_VERSION,
        preprocessing_version = HWSD_CN_PREPROCESSING_VERSION,
        source_version = "HWSD v2.01",
        deep_rule = :extend_deepest_density,
        spatial_policy = "spherical-area mean over HWSD soil pixels",
        conservation = soil_cn_conservation(accumulator),
    )
    raw_targets = finish_soil_cn(
        accumulator, all_cells(grid); minimum_coverage = 0.99, provenance,
    )
    targets, fallback = apply_fallback(
        raw_targets, grid; max_radius_cells = max_fill_radius_cells, production_active,
    )
    write_soil_cn_targets(output_path, targets)
    NCDataset(output_path, "a") do dataset
        defVar(dataset, "fallback_used", Int8, ("cell",))[:] = Int8.(fallback.fallback_used)
        defVar(dataset, "donor_cell_id", Int32, ("cell",))[:] = fallback.donor_cell_id
        distance = defVar(dataset, "fill_distance_km", Float32, ("cell",))
        distance[:] = fallback.fill_distance_km
        distance.attrib["units"] = "km"
        defVar(dataset, "original_minimum_coverage", Float32, ("cell",))[:] =
            fallback.original_minimum_coverage
        defVar(dataset, "soil_area_m2", Float64, ("cell",))[:] = accumulator.target_area
        defVar(dataset, "production_active", Int8, ("cell",))[:] = Int8.(production_active)
        dataset.attrib["production_mask_source"] = production_mask_source
    end
    report = write_qc(
        qc_path, output_path, grid, accumulator, targets, fallback,
        production_mask_source,
    )
    report["active_unresolved_cells"] == 0 || error(
        "$(report["active_unresolved_cells"]) production-active cells remain unresolved; " *
        "increase MAX_FILL_RADIUS_CELLS after reviewing fallback distance QC",
    )
    maximum(report["carbon_relative_error"]) <= 1e-12 || error("carbon conservation QC failed")
    maximum(report["nitrogen_relative_error"]) <= 1e-12 || error("nitrogen conservation QC failed")
    return report
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (4, 5, 6) || error(
        "usage: prepare_canonical_hwsd.jl HWSD_DIR GRID_NC OUTPUT_NC QC_TOML " *
        "[MAX_FILL_RADIUS_CELLS [LANDUSE_NC]]",
    )
    report = prepare_canonical_hwsd(
        abspath(ARGS[1]), abspath(ARGS[2]), abspath(ARGS[3]), abspath(ARGS[4]);
        max_fill_radius_cells = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 8,
        landuse_path = length(ARGS) == 6 ? abspath(ARGS[6]) : nothing,
    )
    println("HWSD canonical-grid preprocessing complete: ", report)
end
