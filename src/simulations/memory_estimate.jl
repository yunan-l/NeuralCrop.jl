function _array_storage_bytes(value, seen::IdDict{Any, Nothing})
    value === nothing && return 0
    if value isa AbstractArray
        haskey(seen, value) && return 0
        seen[value] = nothing
        return length(value) * sizeof(eltype(value))
    end
    value isa Union{Number, AbstractString, Symbol, Function, Type, Module} && return 0
    if value isa Tuple || isstructtype(typeof(value))
        bytes = 0
        for index in 1:fieldcount(typeof(value))
            bytes += _array_storage_bytes(getfield(value, index), seen)
        end
        return bytes
    end
    return 0
end

_array_storage_bytes(value) = _array_storage_bytes(value, IdDict{Any, Nothing}())

function _projected_output_bytes(cells::Int, days::Int, ::Type{T}) where {T}
    annual_rows = fld(days, 365)
    daily_float_fields = length(_DAILY_CROP_FLOAT_OUTPUT_FIELDS) +
        length(_DAILY_SOIL_FLOAT_OUTPUT_FIELDS)
    daily_int_fields = length(_DAILY_CROP_INTEGER_OUTPUT_FIELDS) +
        length(_DAILY_CALENDAR_INTEGER_OUTPUT_FIELDS)
    annual_float_fields = length(_ANNUAL_CROP_FLOAT_OUTPUT_FIELDS)
    annual_int_fields = length(_ANNUAL_CALENDAR_INTEGER_OUTPUT_FIELDS)
    accumulator_bytes = cells * (12 * sizeof(T) + sizeof(Int32))
    return cells * (
        days * (daily_float_fields * sizeof(T) + daily_int_fields * sizeof(Int32)) +
        annual_rows * (annual_float_fields * sizeof(T) + annual_int_fields * sizeof(Int32))
    ) + accumulator_bytes
end

function _stream_output_bytes(
    cells::Int,
    block_days::Int,
    ::Type{T},
    stream::OutputStream,
) where {T}
    daily_bytes = 0
    annual_bytes = 0
    accumulator_bytes = 0
    annual_rows = cld(block_days, 365)
    for variable in stream.variables
        bytes = variable.field in (
            _DAILY_CROP_INTEGER_OUTPUT_FIELDS...,
            _DAILY_CALENDAR_INTEGER_OUTPUT_FIELDS...,
            _ANNUAL_CALENDAR_INTEGER_OUTPUT_FIELDS...,
        ) ? sizeof(Int32) : sizeof(T)
        if _is_annual(variable)
            annual_bytes += annual_rows * cells * bytes
        else
            daily_bytes += block_days * cells * bytes
            stream.frequency !== :daily && (accumulator_bytes += cells * bytes)
        end
    end
    return daily_bytes + annual_bytes + accumulator_bytes
end

# Includes the dynamic-sowing climate calendar: monthly PET and its rolling
# climatology (24 values), annual PET history (365 values), and four Int32
# calendar fields per cell.
const _PERSISTENT_FLOAT_VALUES_PER_CELL = 1636
const _PERSISTENT_NONFLOAT_BYTES_PER_CELL = 46
const _PERSISTENT_FIXED_FLOAT_VALUES = 26

function _estimated_persistent_state_bytes(cells::Int, ::Type{T}) where {T}
    return cells * (
        _PERSISTENT_FLOAT_VALUES_PER_CELL * sizeof(T) +
        _PERSISTENT_NONFLOAT_BYTES_PER_CELL
    ) + _PERSISTENT_FIXED_FLOAT_VALUES * sizeof(T)
end

function _estimated_diagnostics_bytes(cells::Int, days::Int, ::Type{T}) where {T}
    fields = fieldcount(Base.unwrap_unionall(WaterBalance)) +
        fieldcount(Base.unwrap_unionall(NitrogenBalance)) +
        fieldcount(Base.unwrap_unionall(CarbonBalance)) +
        fieldcount(Base.unwrap_unionall(ThermalBalance))
    return fields * cells * days * sizeof(T)
end

function _memory_estimate(
    cells::Int,
    days::Int,
    ::Type{T},
    persistent_state_bytes::Int,
    diagnostics_bytes::Int;
    block_days::Integer,
    prefetch::Bool,
    safety_factor::Real,
    backend::Symbol,
    warmup_years::Integer = 0,
    cached_forcing_blocks::Integer = 0,
    output_stream::Union{Nothing, OutputStream} = nothing,
) where {T}
    block_days > 0 || throw(ArgumentError("block_days must be positive"))
    safety_factor >= 1 || throw(ArgumentError("safety_factor must be at least one"))
    warmup_years >= 0 || throw(ArgumentError("warmup_years must be non-negative"))
    cached_forcing_blocks >= 0 || throw(ArgumentError(
        "cached_forcing_blocks must be non-negative",
    ))
    backend in (:cpu, :accelerator) || throw(ArgumentError(
        "backend must be :cpu or :accelerator",
    ))

    output_bytes = isnothing(output_stream) ?
        _projected_output_bytes(cells, days, T) :
        _stream_output_bytes(cells, Int(block_days), T, output_stream)
    output_growth_bytes = isnothing(output_stream) ?
        days * cells * max(sizeof(T), sizeof(Int32)) : 0
    host_output_chunk_bytes = isnothing(output_stream) ? 0 : output_bytes
    # Six time×cell forcing fields (including nitrate and ammonium deposition)
    # plus one daily global CO₂ vector.
    forcing_block_bytes = (6 * block_days * cells + block_days) * sizeof(T)
    host_forcing_bytes = forcing_block_bytes * (2 + Int(prefetch))
    # Ten annual soil summaries plus carbon/nitrogen target corrections.
    warmup_history_bytes = warmup_years * cells * 12 * sizeof(T)
    cached_forcing_bytes = cached_forcing_blocks * forcing_block_bytes
    model_payload_bytes = persistent_state_bytes + diagnostics_bytes + output_bytes
    backend_peak_bytes = model_payload_bytes + output_growth_bytes + forcing_block_bytes
    host_peak_bytes = backend === :cpu ?
        model_payload_bytes + output_growth_bytes + host_forcing_bytes +
            host_output_chunk_bytes + warmup_history_bytes + cached_forcing_bytes :
        host_forcing_bytes + host_output_chunk_bytes + warmup_history_bytes +
            cached_forcing_bytes
    device_peak_bytes = backend === :cpu ? 0 : backend_peak_bytes

    recommended(value) = ceil(Int, value * safety_factor)
    gib(value) = value / 2.0^30
    return (
        cells,
        days,
        block_days = Int(block_days),
        prefetch,
        backend,
        persistent_state_bytes,
        diagnostics_bytes,
        projected_output_bytes = output_bytes,
        output_growth_bytes,
        host_output_chunk_bytes,
        forcing_block_bytes,
        host_forcing_bytes,
        warmup_history_bytes,
        cached_forcing_bytes,
        model_payload_bytes,
        host_peak_bytes,
        device_peak_bytes,
        recommended_host_peak_bytes = recommended(host_peak_bytes),
        recommended_device_peak_bytes = recommended(device_peak_bytes),
        recommended_host_peak_gib = gib(recommended(host_peak_bytes)),
        recommended_device_peak_gib = gib(recommended(device_peak_bytes)),
        safety_factor = float(safety_factor),
        streaming_output = !isnothing(output_stream),
    )
end

"""
    estimate_memory(cells, days; T=Float32, diagnostics=true, block_days,
                    backend=:accelerator, prefetch=false, safety_factor=1.2)

Conservatively estimate memory before allocating a simulation. The persistent
state coefficients are checked against the current model layout by tests.
"""
function estimate_memory(
    cells::Integer,
    days::Integer;
    T::Type{<:AbstractFloat} = Float32,
    diagnostics::Bool = true,
    block_days::Integer,
    backend::Symbol = :accelerator,
    prefetch::Bool = false,
    safety_factor::Real = 1.2,
    warmup_years::Integer = 0,
    cached_forcing_blocks::Integer = 0,
    output_stream::Union{Nothing, OutputStream} = nothing,
)
    cells > 0 || throw(ArgumentError("cells must be positive"))
    days > 0 || throw(ArgumentError("days must be positive"))
    cell_count = Int(cells)
    day_count = Int(days)
    persistent_state_bytes = _estimated_persistent_state_bytes(cell_count, T)
    diagnostics_bytes = diagnostics ?
        _estimated_diagnostics_bytes(cell_count, day_count, T) : 0
    return _memory_estimate(
        cell_count, day_count, T, persistent_state_bytes, diagnostics_bytes;
        block_days, prefetch, safety_factor, backend, warmup_years,
        cached_forcing_blocks, output_stream,
    )
end

"""
    estimate_memory(simulation; block_days, prefetch=false, safety_factor=1.2)

Estimate model payload memory before running a configured simulation. The
estimate includes output, diagnostics, and climate-transfer buffers.
`prefetch=true` reserves one additional host forcing block. Pass the same
`OutputStream` used by the runner to estimate selected reusable output buffers
instead of retained full time series.

The returned values are bytes. `recommended_*_peak_bytes` apply `safety_factor`
to account for allocator, library, and kernel overhead that cannot be inferred
from model arrays.
"""
function estimate_memory(
    simulation::CropSimulation;
    block_days::Integer,
    prefetch::Bool = false,
    safety_factor::Real = 1.2,
    warmup_years::Integer = 0,
    cached_forcing_blocks::Integer = 0,
    output_stream::Union{Nothing, OutputStream} = nothing,
)
    cells = length(simulation.state.inputs.weather.temp)
    days = simulation.config.days
    T = simulation.config.T
    state_without_output = (
        simulation.processes,
        simulation.state.prognostic,
        simulation.state.fluxes,
        simulation.state.auxiliary,
        simulation.state.inputs,
        simulation.state.events,
        simulation.state.workspace,
    )
    persistent_state_bytes = _array_storage_bytes(state_without_output)
    diagnostics_bytes = _array_storage_bytes(simulation.diagnostics)
    backend_is_host = simulation.config.device === identity
    return _memory_estimate(
        cells, days, T, persistent_state_bytes, diagnostics_bytes;
        block_days, prefetch, safety_factor, warmup_years, cached_forcing_blocks,
        backend = backend_is_host ? :cpu : :accelerator,
        output_stream,
    )
end

"""Estimate memory using the block size and cell selection of a climate reader."""
function estimate_memory(simulation::CropSimulation, climate_reader; kwargs...)
    hasproperty(climate_reader, :block_days) || throw(ArgumentError(
        "climate_reader must expose block_days",
    ))
    if hasproperty(climate_reader, :selection)
        reader_cells = length(climate_reader.selection.compact_indices)
        model_cells = length(simulation.state.inputs.weather.temp)
        reader_cells == model_cells || throw(DimensionMismatch(
            "climate reader has $reader_cells cells, simulation has $model_cells",
        ))
    end
    return estimate_memory(simulation; block_days = climate_reader.block_days, kwargs...)
end
