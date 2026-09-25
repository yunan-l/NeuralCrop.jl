push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))

using NeuralCropData

const SOURCE_RESOLUTION = 1 / 120
const SOURCE_COLUMNS = 43_200
const SOURCE_ROW_BYTES = 2 * SOURCE_COLUMNS
const SOURCE_NODATA = typemax(UInt16)
const EARTH_RADIUS_M = 6_371_000.0

function read_hwsd_tile(path, center_longitude, center_latitude)
    left = center_longitude - 0.25
    bottom = center_latitude - 0.25
    longitude = collect(range(
        left + SOURCE_RESOLUTION / 2; step = SOURCE_RESOLUTION, length = 60,
    ))
    latitude = collect(range(
        bottom + SOURCE_RESOLUTION / 2; step = SOURCE_RESOLUTION, length = 60,
    ))
    mapping_units = Vector{UInt16}(undef, length(longitude) * length(latitude))
    position = 0
    open(path, "r") do io
        for lat in latitude, lon in longitude
            position += 1
            column = floor(Int, (lon + 180) / SOURCE_RESOLUTION) + 1
            row = floor(Int, (90 - lat) / SOURCE_RESOLUTION) + 1
            seek(io, (row - 1) * SOURCE_ROW_BYTES + 2 * (column - 1))
            mapping_units[position] = ltoh(read(io, UInt16))
        end
    end
    latitude_edges = range(
        bottom; step = SOURCE_RESOLUTION, length = length(latitude) + 1,
    )
    longitude_width = deg2rad(SOURCE_RESOLUTION)
    pixel_area = Float64[]
    for latitude_index in eachindex(latitude), _ in eachindex(longitude)
        south = latitude_edges[latitude_index]
        north = latitude_edges[latitude_index + 1]
        push!(
            pixel_area,
            EARTH_RADIUS_M^2 * longitude_width *
                (sin(deg2rad(north)) - sin(deg2rad(south))),
        )
    end
    return (; mapping_units, pixel_area)
end

function query_layer_rows(database_path, mapping_units)
    ids = sort!(Int.(filter(!=(SOURCE_NODATA), unique(mapping_units))))
    isempty(ids) && return Vector{Vector{String}}()
    predicate = join(("HWSD2_SMU_ID=$id" for id in ids), " OR ")
    sql = """
        SELECT HWSD2_SMU_ID,SEQUENCE,SHARE,ROOT_DEPTH,WRB4,LAYER,COARSE,BULK,ORG_CARBON,TOTAL_N
        FROM HWSD2_LAYERS
        WHERE $predicate;
        """
    command = pipeline(
        `mdb-sql -P -H -F -d , $database_path`,
        stdin = IOBuffer(sql),
    )
    lines = filter(!isempty, strip.(readlines(command)))
    isempty(lines) && error("HWSD layer query returned no rows")
    return [split(line, ','; keepempty = true) for line in lines]
end

missing_value(value) = isempty(value) || parse(Float64, value) < 0
parse_value(value) = missing_value(value) ? NaN : parse(Float64, value)

function mapping_unit_targets(rows)
    by_mapping_unit = Dict{Int, Vector{Vector{String}}}()
    for row in rows
        push!(get!(by_mapping_unit, parse(Int, row[1]), Vector{Vector{String}}()), row)
    end
    targets = Dict{Int, NamedTuple}()
    for (mapping_unit, records) in by_mapping_unit
        soil_codes = unique(row[5] for row in records if !isempty(row[5]))
        is_soil = !isempty(soil_codes) && any(code -> code ∉ ("WR", "GG"), soil_codes)
        sequences = sort!(unique(parse(Int, row[2]) for row in records))
        sequence_index = Dict(sequence => index for (index, sequence) in pairs(sequences))
        attributes = [fill(NaN, 7, length(sequences)) for _ in 1:4]
        shares = zeros(Float64, length(sequences))
        root_depth = zeros(Int, length(sequences))
        for row in records
            component = sequence_index[parse(Int, row[2])]
            layer = parse(Int, string(last(row[6])))
            shares[component] = parse(Float64, row[3])
            isempty(row[4]) || (root_depth[component] = parse(Int, row[4]))
            attributes[1][layer, component] = parse_value(row[9]) # organic carbon
            attributes[2][layer, component] = parse_value(row[10]) # total nitrogen
            attributes[3][layer, component] = parse_value(row[8]) # bulk density
            attributes[4][layer, component] = parse_value(row[7]) # coarse fragments
        end
        stocks = hwsd_layer_stocks(attributes...; T = Float64)
        mixed = mix_hwsd_components(
            stocks.carbon, stocks.nitrogen, shares, root_depth; T = Float64,
        )
        remapped_carbon = remap_hwsd_layers(reshape(mixed.carbon, :, 1); T = Float64)
        remapped_nitrogen = remap_hwsd_layers(reshape(mixed.nitrogen, :, 1); T = Float64)
        remapped_uncertain = copy(vec(remapped_carbon.uncertain .| remapped_nitrogen.uncertain))
        for target in axes(NEURALCROP_SOIL_LAYER_BOUNDS, 1),
            source in axes(HWSD_LAYER_BOUNDS, 1)
            overlap = min(
                NEURALCROP_SOIL_LAYER_BOUNDS[target, 2],
                HWSD_LAYER_BOUNDS[source, 2],
            ) - max(
                NEURALCROP_SOIL_LAYER_BOUNDS[target, 1],
                HWSD_LAYER_BOUNDS[source, 1],
            )
            overlap > 0 && mixed.uncertain[source] && (remapped_uncertain[target] = true)
        end
        targets[mapping_unit] = (
            carbon = vec(remapped_carbon.values),
            nitrogen = vec(remapped_nitrogen.values),
            uncertain = remapped_uncertain,
            is_soil,
        )
    end
    return targets
end

function aggregate_hwsd_cell(
    data_directory,
    center_longitude,
    center_latitude,
    cell_id,
)
    raster_path = joinpath(data_directory, "HWSD2.bil")
    database_path = joinpath(data_directory, "HWSD2.mdb")
    tile = read_hwsd_tile(raster_path, center_longitude, center_latitude)
    targets_by_unit = mapping_unit_targets(
        query_layer_rows(database_path, tile.mapping_units),
    )
    pixels = length(tile.mapping_units)
    carbon = fill(NaN, 5, pixels)
    nitrogen = fill(NaN, 5, pixels)
    uncertain = falses(5, pixels)
    for pixel in 1:pixels
        mapping_unit = Int(tile.mapping_units[pixel])
        haskey(targets_by_unit, mapping_unit) || continue
        target = targets_by_unit[mapping_unit]
        carbon[:, pixel] = target.carbon
        nitrogen[:, pixel] = target.nitrogen
        uncertain[:, pixel] = target.uncertain
    end
    aggregator = init_soil_cn_aggregator(Float64, 5, 1)
    target_indices = ifelse.(tile.mapping_units .== SOURCE_NODATA, Int32(0), Int32(1))
    for pixel in 1:pixels
        mapping_unit = Int(tile.mapping_units[pixel])
        target = get(targets_by_unit, mapping_unit, nothing)
        !isnothing(target) && !target.is_soil && (target_indices[pixel] = 0)
    end
    accumulate_soil_cn!(
        aggregator,
        carbon,
        nitrogen,
        target_indices,
        tile.pixel_area;
        uncertain,
    )
    provenance = (
        schema_version = DATA_SCHEMA_VERSION,
        preprocessing_version = HWSD_CN_PREPROCESSING_VERSION,
        source_version = "HWSD v2.01",
        deep_rule = :extend_deepest_density,
        center_longitude = Float64(center_longitude),
        center_latitude = Float64(center_latitude),
        component_policy = "SHARE-weighted; below ROOT_DEPTH is zero; >=99% resolved land SHARE",
        spatial_policy = "area mean over HWSD soil pixels; NODATA, water, and glacier excluded",
        conservation = soil_cn_conservation(aggregator),
    )
    selection = CellSelection([1], [cell_id])
    targets = finish_soil_cn(
        aggregator,
        selection;
        minimum_coverage = 0.99,
        provenance,
    )
    return targets
end

complete_profile(targets) = all(isfinite, targets.soil_organic_carbon) &&
    all(isfinite, targets.total_nitrogen)

function distance_km(longitude1, latitude1, longitude2, latitude2)
    delta_longitude = deg2rad(longitude2 - longitude1)
    delta_latitude = deg2rad(latitude2 - latitude1)
    a = sin(delta_latitude / 2)^2 + cos(deg2rad(latitude1)) *
        cos(deg2rad(latitude2)) * sin(delta_longitude / 2)^2
    return 2 * EARTH_RADIUS_M / 1000 * asin(sqrt(a))
end

function nearest_complete_profile(
    data_directory,
    center_longitude,
    center_latitude;
    cell_id,
    max_radius_cells,
)
    max_radius_cells >= 1 || throw(ArgumentError("max_radius_cells must be positive"))
    candidates = Tuple{Float64, Float64, Float64}[]
    for longitude_offset in -max_radius_cells:max_radius_cells,
        latitude_offset in -max_radius_cells:max_radius_cells
        max(abs(longitude_offset), abs(latitude_offset)) == 0 && continue
        longitude = mod(center_longitude + 0.5 * longitude_offset + 180, 360) - 180
        latitude = center_latitude + 0.5 * latitude_offset
        -89.75 <= latitude <= 89.75 || continue
        distance = distance_km(center_longitude, center_latitude, longitude, latitude)
        push!(candidates, (distance, longitude, latitude))
    end
    sort!(candidates; by = first)
    for (distance, longitude, latitude) in candidates
        donor = aggregate_hwsd_cell(
            data_directory, longitude, latitude, cell_id,
        )
        complete_profile(donor) || continue
        return (; donor, longitude, latitude, distance)
    end
    return nothing
end

function extract_hwsd_cell(
    data_directory,
    center_longitude,
    center_latitude,
    output_path;
    cell_id = 0,
    max_fill_radius_cells = 8,
)
    targets = aggregate_hwsd_cell(
        data_directory, center_longitude, center_latitude, cell_id,
    )
    if !complete_profile(targets)
        original_minimum_coverage = minimum(targets.coverage)
        fill = nearest_complete_profile(
            data_directory, center_longitude, center_latitude;
            cell_id, max_radius_cells = max_fill_radius_cells,
        )
        isnothing(fill) && error(
            "no complete HWSD profile found within $(0.5 * max_fill_radius_cells)° of " *
            "($center_longitude, $center_latitude)",
        )
        provenance = merge(targets.provenance, (
            fill_policy = "nearest complete 0.5-degree HWSD cell",
            donor_longitude = fill.longitude,
            donor_latitude = fill.latitude,
            fill_distance_km = fill.distance,
            original_minimum_coverage = original_minimum_coverage,
        ))
        targets = SoilCNTargets(
            targets.selection,
            fill.donor.layer_bounds,
            fill.donor.soil_organic_carbon,
            fill.donor.total_nitrogen,
            fill.donor.coverage,
            trues(size(fill.donor.uncertain)),
            provenance,
        )
    end
    write_soil_cn_targets(output_path, targets)
    return targets
end

function extract_hwsd_grid_cell(
    data_directory,
    grid_path,
    longitude_index,
    latitude_index,
    output_path;
    cell_id = 0,
    max_fill_radius_cells = 8,
)
    grid = read_grid(grid_path)
    checkbounds(grid.longitude, longitude_index)
    checkbounds(grid.latitude, latitude_index)
    return extract_hwsd_cell(
        data_directory,
        grid.longitude[longitude_index],
        grid.latitude[latitude_index],
        output_path;
        cell_id,
        max_fill_radius_cells,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (5, 6, 7) || error(
        "usage: extract_hwsd_cell.jl DATA_DIR GRID_NC LON_INDEX LAT_INDEX OUTPUT " *
        "[CELL_ID] [MAX_FILL_RADIUS_CELLS]",
    )
    targets = extract_hwsd_grid_cell(
        ARGS[1], ARGS[2], parse(Int, ARGS[3]), parse(Int, ARGS[4]), ARGS[5];
        cell_id = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 0,
        max_fill_radius_cells = length(ARGS) == 7 ? parse(Int, ARGS[7]) : 8,
    )
    println("layer_bounds_m = ", targets.layer_bounds)
    println("soil_organic_carbon_gC_m2 = ", vec(targets.soil_organic_carbon))
    println("total_nitrogen_gN_m2 = ", vec(targets.total_nitrogen))
    println("coverage = ", vec(targets.coverage))
    println("uncertain = ", vec(targets.uncertain))
end
