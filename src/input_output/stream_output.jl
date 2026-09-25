const _STREAM_DAILY_FIELDS = (
    crop = (_DAILY_CROP_FLOAT_OUTPUT_FIELDS..., _DAILY_CROP_INTEGER_OUTPUT_FIELDS...),
    soil = _DAILY_SOIL_FLOAT_OUTPUT_FIELDS,
    calendar = _DAILY_CALENDAR_INTEGER_OUTPUT_FIELDS,
)
const _STREAM_ANNUAL_FIELDS = (
    crop = _ANNUAL_CROP_FLOAT_OUTPUT_FIELDS,
    calendar = _ANNUAL_CALENDAR_INTEGER_OUTPUT_FIELDS,
)
const _MONTH_END_DAYS = cumsum((31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))

struct OutputVariable
    group::Symbol
    field::Symbol
    reduction::Symbol
end

function OutputVariable(group::Symbol, field::Symbol; reduction::Union{Nothing, Symbol} = nothing)
    daily = group in keys(_STREAM_DAILY_FIELDS) && field in getproperty(_STREAM_DAILY_FIELDS, group)
    annual = group in keys(_STREAM_ANNUAL_FIELDS) && field in getproperty(_STREAM_ANNUAL_FIELDS, group)
    (daily || annual) || throw(ArgumentError("unsupported output variable $group.$field"))
    output_variable_spec(group, field)
    default_reduction = field in (
        :gpp, :npp, :respiration, :ecosystem_respiration,
        :heterotrophic_respiration, :evapotranspiration, :sowing_event, :harvest_event,
    ) ?
        :sum : field in (:biomass, :storage_carbon, :growing_mask, :harvesting_mask) ?
        :last : :mean
    selected_reduction = something(reduction, default_reduction)
    selected_reduction in (:sum, :mean, :last) || throw(ArgumentError(
        "reduction must be :sum, :mean, or :last",
    ))
    return OutputVariable(group, field, selected_reduction)
end

struct OutputChunk
    index::Int
    frequency::Symbol
    time::Vector{Int}
    cell_ids::Vector{Int32}
    values::Dict{Symbol, Any}
end

mutable struct OutputStream{W}
    variables::Vector{OutputVariable}
    frequency::Symbol
    writer::W
    cell_ids::Vector{Int32}
    first_output_day::Int
    chunk_index::Int
    period_end_day::Int
    period_count::Int
    accumulators::Dict{Symbol, Any}
    pending_annual::Dict{Int, Dict{Symbol, Any}}
    written_paths::Vector{String}
end

"""
    OutputStream(variables; frequency=:daily, writer, cell_ids, first_output_day=1)

Stream selected diagnostics beginning on the first day of a 365-day simulation
year. The model still advances continuously before `first_output_day`.
"""
function OutputStream(
    variables::AbstractVector{OutputVariable};
    frequency::Symbol = :daily,
    writer = (_ -> nothing),
    cell_ids::AbstractVector{<:Integer},
    first_output_day::Integer = 1,
)
    isempty(variables) && throw(ArgumentError("at least one output variable is required"))
    frequency in (:daily, :monthly, :annual) || throw(ArgumentError(
        "frequency must be :daily, :monthly, or :annual",
    ))
    length(unique((variable.group, variable.field) for variable in variables)) == length(variables) ||
        throw(ArgumentError("output variables must be unique"))
    first_output_day >= 1 || throw(ArgumentError("first_output_day must be positive"))
    (first_output_day - 1) % 365 == 0 || throw(ArgumentError(
        "first_output_day must be the first day of a 365-day simulation year",
    ))
    return OutputStream(
        collect(variables), frequency, writer, Int32.(cell_ids), Int(first_output_day), 0, 0, 0,
        Dict{Symbol, Any}(), Dict{Int, Dict{Symbol, Any}}(), String[],
    )
end

struct JLD2BlockWriter
    directory::String
    prefix::String
end

function JLD2BlockWriter(directory::AbstractString; prefix::AbstractString = "output")
    mkpath(directory)
    return JLD2BlockWriter(String(directory), String(prefix))
end

struct NetCDFBlockWriter
    directory::String
    prefix::String
end

function NetCDFBlockWriter(directory::AbstractString; prefix::AbstractString = "output")
    mkpath(directory)
    return NetCDFBlockWriter(String(directory), String(prefix))
end

function (writer::JLD2BlockWriter)(chunk::OutputChunk)
    path = joinpath(
        writer.directory,
        "$(writer.prefix)_$(chunk.frequency)_$(lpad(chunk.index, 6, '0')).jld2",
    )
    jldsave(
        path;
        frequency = chunk.frequency,
        time = chunk.time,
        cell_ids = chunk.cell_ids,
        values = chunk.values,
    )
    return path
end

function (writer::NetCDFBlockWriter)(chunk::OutputChunk)
    path = joinpath(
        writer.directory,
        "$(writer.prefix)_$(chunk.frequency)_$(lpad(chunk.index, 6, '0')).nc",
    )
    NCDataset(path, "c") do dataset
        defDim(dataset, "time", length(chunk.time))
        defDim(dataset, "cell", length(chunk.cell_ids))
        defVar(dataset, "time", Int32, ("time",))[:] = Int32.(chunk.time)
        defVar(dataset, "cell_id", Int32, ("cell",))[:] = chunk.cell_ids
        for (name, values) in chunk.values
            variable = defVar(dataset, String(name), eltype(values), ("time", "cell"))
            parts = split(String(name), '_'; limit = 2)
            if length(parts) == 2
                spec, _ = output_variable_spec(Symbol(parts[1]), Symbol(parts[2]))
                variable.attrib["units"] = spec.units
                variable.attrib["long_name"] = spec.description
            end
            variable[:, :] = values
        end
    end
    return path
end

_output_key(variable::OutputVariable) = Symbol(variable.group, :_, variable.field)
_is_annual(variable::OutputVariable) =
    variable.group in keys(_STREAM_ANNUAL_FIELDS) &&
    variable.field in getproperty(_STREAM_ANNUAL_FIELDS, variable.group)

function _emit_chunk!(stream::OutputStream, frequency::Symbol, time::Int, values::Dict{Symbol, Any})
    isempty(values) && return nothing
    stream.chunk_index += 1
    chunk = OutputChunk(
        stream.chunk_index, frequency, [time], copy(stream.cell_ids), values,
    )
    result = stream.writer(chunk)
    result isa AbstractString && push!(stream.written_paths, String(result))
    return chunk
end

function _emit_daily_chunk!(
    stream::OutputStream,
    first_day::Int,
    rows::Int,
    values::Dict{Symbol, Any},
)
    isempty(values) && return nothing
    stream.chunk_index += 1
    chunk = OutputChunk(
        stream.chunk_index, :daily, collect(first_day:(first_day + rows - 1)),
        copy(stream.cell_ids), values,
    )
    result = stream.writer(chunk)
    result isa AbstractString && push!(stream.written_paths, String(result))
    return chunk
end

function _period_end(day::Int, frequency::Symbol)
    year_start = fld(day - 1, 365) * 365
    frequency === :annual && return year_start + 365
    day_of_year = day - year_start
    return year_start + _MONTH_END_DAYS[findfirst(>=(day_of_year), _MONTH_END_DAYS)]
end

function _aggregate_row!(
    stream::OutputStream,
    variable::OutputVariable,
    values::AbstractMatrix,
    row::Int,
)
    key = _output_key(variable)
    source = @view values[row, :]
    if !haskey(stream.accumulators, key)
        stream.accumulators[key] = copy(source)
    elseif variable.reduction === :last
        copyto!(stream.accumulators[key], source)
    else
        stream.accumulators[key] .+= source
    end
    return nothing
end

function _flush_period!(stream::OutputStream, end_day::Int)
    values = Dict{Symbol, Any}()
    for variable in stream.variables
        _is_annual(variable) && continue
        key = _output_key(variable)
        haskey(stream.accumulators, key) || continue
        value = stream.accumulators[key]
        reduced = variable.reduction === :mean ? value ./ stream.period_count : value
        values[key] = reshape(copy(reduced), 1, :)
    end
    if stream.frequency === :annual && haskey(stream.pending_annual, end_day)
        merge!(values, pop!(stream.pending_annual, end_day))
    end
    _emit_chunk!(stream, stream.frequency, end_day, values)
    empty!(stream.accumulators)
    stream.period_count = 0
    stream.period_end_day = 0
    return nothing
end

function _host_selected_values(stream::OutputStream, output::Output)
    daily = Dict{Symbol, Any}()
    annual = Dict{Symbol, Any}()
    for variable in stream.variables
        destination = _is_annual(variable) ? annual : daily
        container = getproperty(output, variable.group)
        destination[_output_key(variable)] = Array(getproperty(container, variable.field))
    end
    return daily, annual
end

"""Copy selected output to the host and pass completed chunks to the configured writer."""
function consume_output!(
    stream::OutputStream,
    output::Output,
    first_day::Integer;
    rows::Union{Nothing, Integer} = nothing,
)
    daily, annual = _host_selected_values(stream, output)
    output_rows = isnothing(rows) ?
        (isempty(daily) ? size(output.crop.gpp, 1) : size(first(values(daily)), 1)) : Int(rows)
    length(stream.cell_ids) == size(output.crop.gpp, 2) || throw(DimensionMismatch(
        "stream has $(length(stream.cell_ids)) cell ids, output has $(size(output.crop.gpp, 2)) cells",
    ))

    annual_rows = isempty(annual) ? 0 : size(first(values(annual)), 1)
    annual_end_days = filter(
        day -> day % 365 == 0,
        Int(first_day):(Int(first_day) + output_rows - 1),
    )
    annual_rows == length(annual_end_days) || throw(DimensionMismatch(
        "found $annual_rows annual rows for $(length(annual_end_days)) year boundaries",
    ))
    first_row = max(1, stream.first_output_day - Int(first_day) + 1)
    first_row > output_rows && return stream
    for (row, day) in enumerate(annual_end_days)
        if day >= stream.first_output_day
            stream.pending_annual[day] = Dict(
                key => value[row:row, :] for (key, value) in annual
            )
        end
    end

    if stream.frequency === :daily
        selected_daily = Dict{Symbol, Any}(
            key => value[first_row:output_rows, :] for (key, value) in daily
        )
        _emit_daily_chunk!(stream, Int(first_day) + first_row - 1,
                           output_rows - first_row + 1, selected_daily)
        for day in annual_end_days
            day >= stream.first_output_day || continue
            _emit_chunk!(stream, :annual, day, pop!(stream.pending_annual, day))
        end
        return stream
    end

    for row in first_row:output_rows
        day = Int(first_day) + row - 1
        period_end = _period_end(day, stream.frequency)
        if stream.period_end_day != 0 && stream.period_end_day != period_end
            _flush_period!(stream, stream.period_end_day)
        end
        stream.period_end_day = period_end
        for variable in stream.variables
            _is_annual(variable) || _aggregate_row!(stream, variable, daily[_output_key(variable)], row)
        end
        stream.period_count += 1
        day == period_end && _flush_period!(stream, period_end)
    end

    if stream.frequency === :monthly
        for day in annual_end_days
            day >= stream.first_output_day || continue
            _emit_chunk!(stream, :annual, day, pop!(stream.pending_annual, day))
        end
    end
    return stream
end

"""Emit the final incomplete month or year, if one remains."""
function finish_output_stream!(stream::OutputStream, final_day::Integer)
    stream.period_count > 0 && _flush_period!(stream, Int(final_day))
    for day in sort!(collect(keys(stream.pending_annual)))
        _emit_chunk!(stream, :annual, day, pop!(stream.pending_annual, day))
    end
    return stream
end

"""Release all emitted time-series rows while preserving annual accumulators."""
function clear_output_timeseries!(output::Output)
    for group in (:crop, :soil, :climate, :calendar)
        container = getproperty(output, group)
        for field in fieldnames(typeof(container))
            values = getproperty(container, field)
            setproperty!(container, field, similar(values, 0, size(values, 2)))
        end
    end
    return output
end

function _output_timeseries_empty(output::Output)
    return all(group -> all(
        field -> size(getproperty(getproperty(output, group), field), 1) == 0,
        fieldnames(typeof(getproperty(output, group))),
    ), (:crop, :soil, :climate, :calendar))
end
