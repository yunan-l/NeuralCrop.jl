const _CLIMATE_DATASETS = (:temp, :prec, :lwnet, :swdown)
const _NITROGEN_DEPOSITION_DATASETS = (:no3_deposition, :nh4_deposition)
const _NOLEAP_MONTH_LENGTHS = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
const _LPJML_MONTH_CENTRE_INTERVALS = (30, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

"""Lazy, restartable iterator over bounded daily climate blocks."""
struct ClimateBlockReader{T <: AbstractFloat, TT, G, S, C}
    catalog::C
    grid::G
    selection::S
    time::Vector{TT}
    indices::Vector{Int}
    block_days::Int
    co2::CO2Series{T}
    calendar::String
    nitrogen_deposition_year::Union{Nothing, Int}
end

"""Read a two-column `year ppm` global annual CO₂ text file."""
function read_co2_series(path::AbstractString; T::Type{<:AbstractFloat} = Float32)
    years = Int32[]
    values = T[]
    for (line_number, line) in enumerate(eachline(path))
        content = strip(first(split(line, '#'; limit = 2)))
        isempty(content) && continue
        fields = split(content)
        length(fields) >= 2 || throw(ArgumentError("invalid CO₂ row $line_number in $path"))
        year = tryparse(Int32, fields[1])
        value = tryparse(T, fields[2])
        isnothing(year) && throw(ArgumentError("invalid CO₂ year on row $line_number in $path"))
        isnothing(value) && throw(ArgumentError("invalid CO₂ value on row $line_number in $path"))
        value > zero(T) || throw(ArgumentError("CO₂ values must be positive"))
        push!(years, year)
        push!(values, value)
    end
    isempty(years) && throw(ArgumentError("CO₂ file contains no data"))
    issorted(years) || throw(ArgumentError("CO₂ years must be sorted"))
    allunique(years) || throw(ArgumentError("CO₂ years must be unique"))
    provenance = DataProvenance(DATA_SCHEMA_VERSION, abspath(path), "co2", "ppm")
    return CO2Series(years, values, provenance)
end

function _climate_time(catalog::DatasetCatalog)
    reference = _time_coordinate(dataset(catalog, :temp), Colon())
    metadata = _time_metadata(dataset(catalog, :temp))
    isempty(reference) && throw(ArgumentError("climate time coordinate cannot be empty"))
    for name in _CLIMATE_DATASETS[2:end]
        candidate = _time_coordinate(dataset(catalog, name), Colon())
        candidate == reference || throw(ArgumentError("$name time coordinate does not match temp"))
        _time_metadata(dataset(catalog, name)) == metadata ||
            throw(ArgumentError("$name time metadata do not match temp"))
    end
    return reference, metadata
end

function _climate_indices(time, start_year, end_year)
    if isnothing(start_year) && isnothing(end_year)
        return 1:length(time)
    end
    first_year = isnothing(start_year) ? minimum(_calendar_year.(time)) : Int(start_year)
    last_year = isnothing(end_year) ? maximum(_calendar_year.(time)) : Int(end_year)
    first_year <= last_year || throw(ArgumentError("start_year must not exceed end_year"))
    indices = findall(value -> first_year <= _calendar_year(value) <= last_year, time)
    isempty(indices) && throw(ArgumentError("requested climate years are unavailable"))
    indices == collect(first(indices):last(indices)) ||
        throw(ArgumentError("requested climate rows must be contiguous"))
    return collect(first(indices):last(indices))
end

const _STANDARD_CALENDARS = ("standard", "gregorian", "proleptic_gregorian")
const _NOLEAP_CALENDARS = ("noleap", "365_day")

function _date_month_day(value)
    value isa Real && return nothing
    try
        return Dates.month(value), Dates.day(value)
    catch
        return nothing
    end
end

function _normalize_calendar_indices(time, indices::AbstractVector{<:Integer}, calendar::AbstractString)
    normalized = lowercase(String(calendar))
    normalized in ("360_day", "360") && throw(ArgumentError(
        "360-day climate is unsupported; NeuralCrop currently requires 365-day years",
    ))
    isempty(normalized) || normalized in _STANDARD_CALENDARS || normalized in _NOLEAP_CALENDARS ||
        throw(ArgumentError("unsupported climate calendar '$calendar'"))
    dates_available = all(!isnothing(_date_month_day(time[index])) for index in indices)
    dates_available || return Int.(indices), isempty(normalized) ? "unspecified" : normalized

    kept = if normalized in _STANDARD_CALENDARS || isempty(normalized)
        [index for index in indices if _date_month_day(time[index]) != (2, 29)]
    else
        Int.(indices)
    end
    counts = Dict{Int, Int}()
    for index in kept
        year = _calendar_year(time[index])
        counts[year] = get(counts, year, 0) + 1
    end
    all(==(365), values(counts)) || throw(ArgumentError(
        "daily climate must contain complete 365-day years after calendar normalization",
    ))
    return kept, normalized in _NOLEAP_CALENDARS ? normalized : "noleap"
end

"""Create a lazy iterator that reads at most `block_days` daily rows at once."""
function climate_blocks(
    catalog::DatasetCatalog,
    grid::GridIndex;
    selection::CellSelection = all_cells(grid),
    co2_path::AbstractString = dataset(catalog, :co2).path,
    start_year = nothing,
    end_year = nothing,
    block_days::Integer = 31,
    T::Type{<:AbstractFloat} = Float32,
    nitrogen_deposition_year::Union{Nothing, Integer} = nothing,
)
    block_days > 0 || throw(ArgumentError("block_days must be positive"))
    time, time_metadata = _climate_time(catalog)
    requested_indices = _climate_indices(time, start_year, end_year)
    indices, calendar = _normalize_calendar_indices(
        time, requested_indices, time_metadata.calendar,
    )
    co2 = read_co2_series(co2_path; T)
    requested_years = unique(Int32.(_calendar_year.(time[indices])))
    available_years = Set(co2.years)
    all(year -> year in available_years, requested_years) ||
        throw(ArgumentError("CO₂ data do not cover all requested climate years"))
    return ClimateBlockReader{T, eltype(time), typeof(grid), typeof(selection), typeof(catalog)}(
        catalog, grid, selection, time, indices, Int(block_days), co2, calendar,
        isnothing(nitrogen_deposition_year) ? nothing : Int(nitrogen_deposition_year),
    )
end

Base.length(reader::ClimateBlockReader) = cld(length(reader.indices), reader.block_days)
climate_days(reader::ClimateBlockReader) = length(reader.indices)

function _block_indices(reader::ClimateBlockReader, block_index::Integer)
    1 <= block_index <= length(reader) || throw(BoundsError(reader, block_index))
    first_position = (block_index - 1) * reader.block_days + 1
    last_position = min(length(reader.indices), first_position + reader.block_days - 1)
    return reader.indices[first_position:last_position]
end

function _validate_climate(name::Symbol, values)
    any(ismissing, values) && throw(ArgumentError("$name contains missing values"))
    all(isfinite, values) || throw(ArgumentError("$name contains non-finite values"))
    name === :prec && !all(>=(0), values) && throw(ArgumentError("precipitation must be non-negative"))
    name === :swdown && !all(>=(0), values) && throw(ArgumentError("shortwave radiation must be non-negative"))
    name === :temp && !all(value -> -100 <= value <= 70, values) &&
        throw(ArgumentError("temperature must be supplied in degrees Celsius"))
    return values
end

function _normalize_climate_units(name::Symbol, values::Matrix{T}, units::AbstractString) where {T}
    normalized = lowercase(replace(strip(units), " " => ""))
    if name === :temp
        if normalized in ("k", "kelvin")
            return values .- T(273.15), "degC"
        elseif normalized in (
            "degc", "°c", "c", "celsius", "degree_celsius", "degrees_celsius",
        )
            return values, "degC"
        end
    elseif name === :prec
        if normalized in ("kgm-2s-1", "kgm^-2s^-1", "kg/m2/s", "mms-1", "mm/s")
            return values .* T(86400), "mm/day"
        elseif normalized in (
            "mm/day", "mmd-1", "mmday-1", "mm",
            "kg/m2/day", "kgm-2day-1", "kgm^-2day^-1",
        )
            return values, "mm/day"
        end
    elseif name in (:lwnet, :swdown)
        if normalized in ("w/m2", "wm-2", "wm^-2")
            return values, "W/m2"
        end
    elseif name in _NITROGEN_DEPOSITION_DATASETS
        if normalized in (
            "g/m2/day", "gm-2day-1", "gm^-2day^-1", "gm-2d-1", "gm^-2d^-1",
        )
            return values, "g/m2/day"
        end
    end
    throw(ArgumentError("unsupported $name units '$units'"))
end

function _monthly_deposition_positions(spec::DatasetSpec, years::Vector{Int})
    source_time = _time_coordinate(spec, Colon())
    all(!isnothing(_date_month_day(value)) for value in source_time) || throw(ArgumentError(
        "monthly nitrogen deposition requires calendar dates in $(spec.path)",
    ))
    positions = Int[]
    for year in years
        year_positions = findall(value -> _calendar_year(value) == year, source_time)
        length(year_positions) == 12 || throw(ArgumentError(
            "monthly nitrogen deposition must provide 12 values for year $year in $(spec.path)",
        ))
        months = first.(_date_month_day.(source_time[year_positions]))
        months == collect(1:12) || throw(ArgumentError(
            "monthly nitrogen deposition must be ordered from January through December in $(spec.path)",
        ))
        append!(positions, year_positions)
    end
    return positions
end

"""Interpolate monthly N deposition to daily values with LPJmL's month-centred rule."""
function _daily_nitrogen_deposition(
    reader::ClimateBlockReader{T}, name::Symbol, block_time,
) where {T}
    haskey(reader.catalog.datasets, name) || return nothing
    all(!isnothing(_date_month_day(value)) for value in block_time) || throw(ArgumentError(
        "daily climate dates are required when nitrogen deposition is configured",
    ))
    years = isnothing(reader.nitrogen_deposition_year) ?
        unique(Int.(_calendar_year.(block_time))) : [reader.nitrogen_deposition_year]
    spec = dataset(reader.catalog, name)
    positions = _monthly_deposition_positions(spec, years)
    monthly = read_compact_variable(
        spec, reader.grid;
        selection = reader.selection,
        selectors = (time = positions,),
        order = (:time, :cell),
        T,
    )
    monthly_values, _ = _normalize_climate_units(
        name, Matrix{T}(monthly.values), monthly.provenance.units,
    )
    year_offsets = Dict(year => 12 * (index - 1) for (index, year) in pairs(years))
    daily = Matrix{T}(undef, length(block_time), size(monthly_values, 2))
    for day_index in eachindex(block_time)
        year = isnothing(reader.nitrogen_deposition_year) ?
            _calendar_year(block_time[day_index]) : reader.nitrogen_deposition_year
        month, day_of_month = _date_month_day(block_time[day_index])
        month_length = _NOLEAP_MONTH_LENGTHS[month]
        previous_month = month == 1 ? 12 : month - 1
        next_month = month == 12 ? 1 : month + 1
        source_month, target_month, offset, denominator = if day_of_month >= month_length ÷ 2
            month, next_month, day_of_month - month_length ÷ 2,
            _LPJML_MONTH_CENTRE_INTERVALS[month]
        else
            previous_month, month,
            day_of_month + (_NOLEAP_MONTH_LENGTHS[previous_month] + 1) ÷ 2,
            _LPJML_MONTH_CENTRE_INTERVALS[previous_month]
        end
        year_offset = year_offsets[year]
        fraction = T(offset) / T(denominator)
        @views daily[day_index, :] .= monthly_values[year_offset + source_month, :] .+
            fraction .* (monthly_values[year_offset + target_month, :] .-
                         monthly_values[year_offset + source_month, :])
    end
    return daily
end

function _daily_co2(series::CO2Series{T}, time) where {T}
    positions = Dict(year => index for (index, year) in pairs(series.years))
    return T[series.values[positions[Int32(_calendar_year(value))]] for value in time]
end

"""Read one numbered block from a `ClimateBlockReader`."""
function read_climate_block(reader::ClimateBlockReader{T}, block_index::Integer) where {T}
    indices = _block_indices(reader, block_index)
    fields = map(_CLIMATE_DATASETS) do name
        read_compact_variable(
            dataset(reader.catalog, name), reader.grid;
            selection = reader.selection,
            selectors = (time = indices,),
            order = (:time, :cell),
            T,
        )
    end
    temp, prec, lwnet, swdown = fields
    temp_values, temp_units = _normalize_climate_units(
        :temp, Matrix{T}(temp.values), temp.provenance.units,
    )
    prec_values, prec_units = _normalize_climate_units(
        :prec, Matrix{T}(prec.values), prec.provenance.units,
    )
    lwnet_values, lwnet_units = _normalize_climate_units(
        :lwnet, Matrix{T}(lwnet.values), lwnet.provenance.units,
    )
    swdown_values, swdown_units = _normalize_climate_units(
        :swdown, Matrix{T}(swdown.values), swdown.provenance.units,
    )
    _validate_climate(:temp, temp_values)
    _validate_climate(:prec, prec_values)
    _validate_climate(:lwnet, lwnet_values)
    _validate_climate(:swdown, swdown_values)
    block_time = collect(reader.time[indices])
    no3_deposition = _daily_nitrogen_deposition(reader, :no3_deposition, block_time)
    nh4_deposition = _daily_nitrogen_deposition(reader, :nh4_deposition, block_time)
    provenance = (
        temp = temp.provenance,
        prec = prec.provenance,
        lwnet = lwnet.provenance,
        swdown = swdown.provenance,
        co2 = reader.co2.provenance,
        no3_deposition = isnothing(no3_deposition) ? nothing : dataset(reader.catalog, :no3_deposition).path,
        nh4_deposition = isnothing(nh4_deposition) ? nothing : dataset(reader.catalog, :nh4_deposition).path,
        model_units = (
            temp = temp_units,
            prec = prec_units,
            lwnet = lwnet_units,
            swdown = swdown_units,
            co2 = "ppm",
            no3_deposition = isnothing(no3_deposition) ? nothing : "g/m2/day",
            nh4_deposition = isnothing(nh4_deposition) ? nothing : "g/m2/day",
        ),
        calendar = reader.calendar,
    )
    return ClimateBlock(
        block_time,
        temp_values,
        prec_values,
        swdown_values,
        lwnet_values,
        no3_deposition,
        nh4_deposition,
        _daily_co2(reader.co2, block_time),
        reader.selection,
        provenance,
    )
end

function Base.iterate(reader::ClimateBlockReader, block_index::Int = 1)
    block_index > length(reader) && return nothing
    return read_climate_block(reader, block_index), block_index + 1
end

"""Return the model-facing forcing tuple without provenance metadata."""
function climate_forcing(block::ClimateBlock)
    forcing = (
        temp = block.temperature,
        prec = block.precipitation,
        sw = block.shortwave,
        lw = block.longwave,
        co2 = block.co2,
        co2_daily = true,
        backend_neutral = true,
    )
    !isnothing(block.no3_deposition) && (forcing = merge(forcing, (
        no3_deposition = block.no3_deposition,
    )))
    !isnothing(block.nh4_deposition) && (forcing = merge(forcing, (
        nh4_deposition = block.nh4_deposition,
    )))
    return forcing
end

"""Expose climate blocks through NeuralCrop's existing block-runner API."""
climate_forcings(source::ClimateBlockReader) = ClimateForcingReader(source)
climate_forcings(source::AbstractVector{<:ClimateBlock}) = ClimateForcingReader(source)

Base.size(forcings::ClimateForcingReader) = (length(forcings.source),)
Base.IndexStyle(::Type{<:ClimateForcingReader}) = IndexLinear()
_climate_block(source::ClimateBlockReader, index::Int) = read_climate_block(source, index)
_climate_block(source::AbstractVector{<:ClimateBlock}, index::Int) = source[index]
Base.getindex(forcings::ClimateForcingReader, index::Int) =
    climate_forcing(_climate_block(forcings.source, index))

"""
    prefetch_climate_forcings(reader)

Return a sequential forcing vector that overlaps model work on block `n` with
reading block `n + 1`. At least two Julia threads are required for true I/O and
compute overlap; the interface remains correct with one thread.
"""
function prefetch_climate_forcings(source::ClimateBlockReader)
    return PrefetchedClimateForcingReader(source, 0, nothing, false)
end

Base.size(forcings::PrefetchedClimateForcingReader) = (length(forcings.source),)
Base.IndexStyle(::Type{<:PrefetchedClimateForcingReader}) = IndexLinear()

function _schedule_prefetch!(forcings::PrefetchedClimateForcingReader, index::Int)
    if index <= length(forcings)
        forcings.scheduled_index = index
        forcings.task = Threads.@spawn climate_forcing(read_climate_block(forcings.source, index))
    else
        forcings.scheduled_index = 0
        forcings.task = nothing
    end
    return forcings
end

function Base.getindex(forcings::PrefetchedClimateForcingReader, index::Int)
    forcings.closed && throw(ArgumentError("prefetched climate reader is closed"))
    checkbounds(forcings, index)
    if isnothing(forcings.task) || forcings.scheduled_index != index
        !isnothing(forcings.task) && fetch(forcings.task)
        _schedule_prefetch!(forcings, index)
    end
    block = fetch(forcings.task)
    _schedule_prefetch!(forcings, index + 1)
    return block
end

function Base.close(forcings::PrefetchedClimateForcingReader)
    if !forcings.closed
        !isnothing(forcings.task) && fetch(forcings.task)
        forcings.task = nothing
        forcings.closed = true
    end
    return nothing
end
