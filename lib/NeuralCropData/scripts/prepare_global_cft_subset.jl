using Dates
using NCDatasets
using TOML

const CFT_DIMENSIONS = Set(("pft", "cft", "crop", "band"))
const TIME_DIMENSIONS = Set(("time", "year"))
dimension_kind(name) = lowercase(String(name)) in CFT_DIMENSIONS ? :cft :
    lowercase(String(name)) in TIME_DIMENSIONS ? :time : :other

attributes(value) = Dict(String(key) => item for (key, item) in pairs(value.attrib)
    if String(key) != "_NCProperties")

function copy_attributes!(destination, source)
    for (key, value) in attributes(source)
        destination.attrib[key] = value
    end
    return destination
end

calendar_year(value::Real) = round(Int, value)
calendar_year(value) = Dates.year(value)

function validate_365_day_years(time, years, label)
    year_values = calendar_year.(time)
    for year in years
        indices = findall(==(year), year_values)
        length(indices) == 365 || throw(ArgumentError(
            "$label year $year contains $(length(indices)) rows; expected exactly 365",
        ))
        if !isempty(indices) && !(time[first(indices)] isa Real)
            any(index -> (Dates.month(time[index]), Dates.day(time[index])) == (2, 29), indices) &&
                throw(ArgumentError("$label year $year contains 29 February"))
            (Dates.month(time[first(indices)]), Dates.day(time[first(indices)])) == (1, 1) ||
                throw(ArgumentError("$label year $year does not start on 1 January"))
            (Dates.month(time[last(indices)]), Dates.day(time[last(indices)])) == (12, 31) ||
                throw(ArgumentError("$label year $year does not end on 31 December"))
        end
    end
    return nothing
end

function indices_for_years(
    dataset,
    variable_name::AbstractString,
    years;
    require_365_days::Bool = true,
)
    dimensions = String.(dimnames(dataset[variable_name]))
    time_position = findfirst(name -> dimension_kind(name) === :time, dimensions)
    isnothing(time_position) && throw(ArgumentError("$variable_name has no time dimension"))
    time_name = dimensions[time_position]
    time = collect(dataset[time_name][:])
    indices = findall(value -> calendar_year(value) in years, time)
    found = unique(calendar_year.(time[indices]))
    found == collect(years) || throw(ArgumentError(
        "$variable_name contains years $found, expected $(collect(years))",
    ))
    require_365_days && validate_365_day_years(time, years, variable_name)
    return time_name, indices
end

function daily_indices_for_years(
    dataset,
    variable_name::AbstractString,
    years,
    source_start_year::Integer,
)
    dimensions = String.(dimnames(dataset[variable_name]))
    time_position = findfirst(name -> dimension_kind(name) === :time, dimensions)
    isnothing(time_position) && throw(ArgumentError("$variable_name has no time dimension"))
    time_name = dimensions[time_position]
    time_length = size(dataset[variable_name], time_position)
    indices = Int[]
    for year in years
        first_index = (year - source_start_year) * 365 + 1
        last_index = first_index + 364
        1 <= first_index <= last_index <= time_length || throw(ArgumentError(
            "$variable_name year $year maps to rows $first_index:$last_index, " *
            "outside its $time_length-row time dimension",
        ))
        append!(indices, first_index:last_index)
    end
    return time_name, indices
end

function selected_time_values(
    path,
    variable_name,
    years;
    daily_source_start_year::Union{Nothing, Integer} = nothing,
)
    return NCDataset(path, "r") do dataset
        time_name, indices = isnothing(daily_source_start_year) ?
            indices_for_years(dataset, variable_name, years) :
            daily_indices_for_years(
                dataset, variable_name, years, daily_source_start_year,
            )
        collect(dataset[time_name][indices])
    end
end

function define_subset_variable!(destination, source, name, selections)
    source_variable = variable(source, name)
    dimensions = String.(dimnames(source_variable))
    output = defVar(
        destination, name, eltype(source_variable), Tuple(dimensions);
        attrib = attributes(source_variable),
    )
    source_indices = ntuple(position -> selections[dimensions[position]], length(dimensions))
    output_indices = ntuple(position -> 1:length(source_indices[position]), length(dimensions))
    output[output_indices...] = source_variable[source_indices...]
    return output
end

function copy_main_variable!(destination, source, name, selections; chunk_length::Integer)
    source_variable = variable(source, name)
    dimensions = String.(dimnames(source_variable))
    output = defVar(
        destination, name, eltype(source_variable), Tuple(dimensions);
        attrib = attributes(source_variable), deflatelevel = 1, shuffle = true,
    )
    time_position = findfirst(dimension -> dimension_kind(dimension) === :time, dimensions)
    if isnothing(time_position)
        source_indices = ntuple(position -> selections[dimensions[position]], length(dimensions))
        output_indices = ntuple(position -> 1:length(source_indices[position]), length(dimensions))
        output[output_indices...] = source_variable[source_indices...]
        return output
    end

    selected_time = selections[dimensions[time_position]]
    for first_output in 1:chunk_length:length(selected_time)
        last_output = min(length(selected_time), first_output + chunk_length - 1)
        output_time = first_output:last_output
        source_indices = ntuple(length(dimensions)) do position
            position == time_position ? selected_time[output_time] : selections[dimensions[position]]
        end
        output_indices = ntuple(length(dimensions)) do position
            position == time_position ? output_time : 1:length(source_indices[position])
        end
        output[output_indices...] = source_variable[source_indices...]
    end
    return output
end

function subset_netcdf(
    input_path::AbstractString,
    output_path::AbstractString,
    variable_name::AbstractString;
    cft_index::Union{Nothing, Integer} = nothing,
    cft_indices::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    years::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    require_365_days::Bool = true,
    daily_source_start_year::Union{Nothing, Integer} = nothing,
    chunk_length::Integer = 31,
)
    chunk_length > 0 || throw(ArgumentError("chunk_length must be positive"))
    isnothing(cft_index) || isnothing(cft_indices) || throw(ArgumentError(
        "specify only one of cft_index or cft_indices",
    ))
    mkpath(dirname(output_path))
    NCDataset(input_path, "r") do source
        haskey(source, variable_name) || throw(ArgumentError(
            "variable '$variable_name' not found in $input_path",
        ))
        dimensions = String.(dimnames(source[variable_name]))
        selections = Dict{String, Vector{Int}}()
        for (position, dimension) in pairs(dimensions)
            selections[dimension] = collect(1:size(source[variable_name], position))
        end
        selected_cft_indices = isnothing(cft_indices) ?
            (isnothing(cft_index) ? nothing : [Int(cft_index)]) : Int.(cft_indices)
        cft_dimension = nothing
        if !isnothing(selected_cft_indices)
            cft_dimensions = filter(dimension -> dimension_kind(dimension) === :cft, dimensions)
            length(cft_dimensions) == 1 || throw(ArgumentError(
                "$variable_name must contain exactly one pft/cft/crop/band dimension",
            ))
            cft_dimension = only(cft_dimensions)
            isempty(selected_cft_indices) && throw(ArgumentError("cft_indices must not be empty"))
            all(index -> 1 <= index <= length(selections[cft_dimension]), selected_cft_indices) ||
                throw(BoundsError(selections[cft_dimension], selected_cft_indices))
            selections[cft_dimension] = selected_cft_indices
        end
        if !isnothing(years)
            time_dimension, time_indices = isnothing(daily_source_start_year) ?
                indices_for_years(source, variable_name, years; require_365_days) :
                daily_indices_for_years(
                    source, variable_name, years, daily_source_start_year,
                )
            selections[time_dimension] = time_indices
        end

        isfile(output_path) && throw(ArgumentError("output already exists: $output_path"))
        NCDataset(output_path, "c") do destination
            copy_attributes!(destination, source)
            destination.attrib["neuralcrop_subset_source"] = abspath(input_path)
            destination.attrib["neuralcrop_subset_variable"] = variable_name
            !isnothing(selected_cft_indices) &&
                (destination.attrib["neuralcrop_source_cft_indices"] = join(selected_cft_indices, ","))
            !isnothing(years) &&
                (destination.attrib["neuralcrop_selected_years"] = join(years, ","))

            for dimension in dimensions
                defDim(destination, dimension, length(selections[dimension]))
            end
            for dimension in dimensions
                haskey(source, dimension) || continue
                coordinate_dimensions = String.(dimnames(source[dimension]))
                all(haskey(selections, name) for name in coordinate_dimensions) || continue
                define_subset_variable!(destination, source, dimension, selections)
                dimension == cft_dimension &&
                    (destination[dimension][:] = collect(1:length(selected_cft_indices)))
            end
            copy_main_variable!(
                destination, source, variable_name, selections; chunk_length,
            )
        end
    end
    return output_path
end

function has_time_dimension(path, variable_name)
    return NCDataset(path, "r") do dataset
        dimensions = String.(dimnames(dataset[variable_name]))
        any(dimension -> dimension_kind(dimension) === :time, dimensions)
    end
end

function management_years(path, variable_name, requested_years)
    lowercase(String(variable_name)) == "sdate" && return nothing
    has_time_dimension(path, variable_name) || return nothing
    isnothing(requested_years) && return nothing
    return requested_years isa Integer ? [Int(requested_years)] : Int.(requested_years)
end

function requested_management_years(settings)
    value = get(settings, "management_years", get(settings, "management_year", 2015))
    value isa AbstractString && lowercase(strip(value)) == "all" && return nothing
    value isa Integer && return [Int(value)]
    value isa AbstractVector && return Int.(value)
    throw(ArgumentError("subset.management_years must be an integer, an array, or 'all'"))
end

function subset_co2(input_path, output_path, years)
    isfile(output_path) && throw(ArgumentError("output already exists: $output_path"))
    mkpath(dirname(output_path))
    selected = String[]
    for line in eachline(input_path)
        content = strip(first(split(line, '#'; limit = 2)))
        isempty(content) && continue
        fields = split(content)
        length(fields) >= 2 || continue
        year = tryparse(Int, fields[1])
        !isnothing(year) && year in years && push!(selected, line)
    end
    found = sort([parse(Int, first(split(line))) for line in selected])
    found == sort(collect(years)) || throw(ArgumentError(
        "CO₂ file contains years $found, expected $(sort(collect(years)))",
    ))
    open(output_path, "w") do output
        println(output, "# Extracted from $(abspath(input_path))")
        foreach(line -> println(output, line), selected)
    end
    return output_path
end

function resolve_path(config_path, value)
    path = String(value)
    return isabspath(path) ? path : normpath(joinpath(dirname(config_path), path))
end

function prepare_subset(config_path::AbstractString)
    config = TOML.parsefile(config_path)
    settings = config["subset"]
    output_directory = resolve_path(config_path, settings["output_directory"])
    management_year_selection = requested_management_years(settings)
    climate_start_year = Int(get(settings, "climate_start_year", 2015))
    climate_year_setting = get(settings, "climate_years", [2015, 2016])
    climate_years = climate_year_setting isa Integer ?
        collect(climate_start_year:(climate_start_year + Int(climate_year_setting) - 1)) :
        Int.(climate_year_setting)
    climate_source_start_year = Int(get(settings, "climate_source_start_year", 1901))
    chunk_length = Int(get(settings, "chunk_length", 31))

    for name in sort!(collect(keys(config["management"])))
        spec = config["management"][name]
        source_cft_indices = Int.(spec["cft_indices"])
        expected_count = name == "residue_fraction" ? 12 : 24
        length(source_cft_indices) == expected_count || throw(ArgumentError(
            "management.$name cft_indices must contain $expected_count entries",
        ))
        input_path = resolve_path(config_path, spec["input"])
        subset_netcdf(
            input_path,
            joinpath(output_directory, spec["output"]),
            spec["variable"];
            cft_indices = source_cft_indices,
            years = management_years(
                input_path, spec["variable"], management_year_selection,
            ),
            require_365_days = false,
            chunk_length,
        )
    end

    if haskey(config, "climate")
        climate_names = sort!(collect(keys(config["climate"])))
        isempty(climate_names) && throw(ArgumentError("climate configuration cannot be empty"))
        reference = config["climate"][first(climate_names)]
        reference_path = resolve_path(config_path, reference["input"])
        selected_years = climate_years
        reference_time = selected_time_values(
            reference_path, reference["variable"], selected_years;
            daily_source_start_year = climate_source_start_year,
        )
        for name in climate_names
            spec = config["climate"][name]
            input_path = resolve_path(config_path, spec["input"])
            selected_time_values(
                input_path, spec["variable"], selected_years;
                daily_source_start_year = climate_source_start_year,
            ) == reference_time || throw(ArgumentError(
                "climate time coordinate for $name does not match the reference",
            ))
            subset_netcdf(
                input_path,
                joinpath(output_directory, spec["output"]),
                spec["variable"];
                years = selected_years,
                daily_source_start_year = climate_source_start_year,
                chunk_length,
            )
        end
        if haskey(config, "co2")
            subset_co2(
                resolve_path(config_path, config["co2"]["input"]),
                joinpath(output_directory, config["co2"]["output"]),
                selected_years,
            )
        end
    elseif haskey(config, "co2")
        throw(ArgumentError("co2 extraction requires a climate year selection"))
    end
    println(
        "Prepared 12 CFTs × rainfed/irrigated management for " *
        (isnothing(management_year_selection) ? "all available years " :
         "years $(join(management_year_selection, ", ")) ") *
        "(12 shared residue bands)" *
        (haskey(config, "climate") ?
         ", and climate years $(join(climate_years, ", "))" : "") *
        " in $output_directory",
    )
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: julia --project=lib/NeuralCropData $(basename(@__FILE__)) CONFIG.toml")
    prepare_subset(abspath(only(ARGS)))
end
