const _NOLEAP_MONTH_LENGTHS = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
const _LPJML_MONTH_CENTRE_INTERVALS = (30, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

"""Expand monthly no-leap forcing using LPJmL's month-centred interpolation."""
function _lpjml_monthly_to_daily(values::AbstractMatrix{T}, daily_days::Integer) where {T <: AbstractFloat}
    daily_days % 365 == 0 || throw(DimensionMismatch(
        "monthly forcing requires a whole number of 365-day years, got $daily_days days",
    ))
    size(values, 1) % 12 == 0 || throw(DimensionMismatch(
        "monthly forcing must have a multiple of 12 rows, got $(size(values, 1))",
    ))
    monthly_years = div(size(values, 1), 12)
    daily_years = div(daily_days, 365)
    monthly_years == daily_years || throw(DimensionMismatch(
        "monthly forcing has $monthly_years year(s), but daily climate has $daily_years",
    ))

    daily = Matrix{T}(undef, daily_days, size(values, 2))
    day = 1
    for year in 0:(daily_years - 1), month in 1:12
        month_length = _NOLEAP_MONTH_LENGTHS[month]
        previous_month = month == 1 ? 12 : month - 1
        next_month = month == 12 ? 1 : month + 1
        for day_of_month in 1:month_length
            source_month, target_month, offset, denominator = if day_of_month >= month_length ÷ 2
                month, next_month, day_of_month - month_length ÷ 2,
                _LPJML_MONTH_CENTRE_INTERVALS[month]
            else
                previous_month, month,
                day_of_month + (_NOLEAP_MONTH_LENGTHS[previous_month] + 1) ÷ 2,
                _LPJML_MONTH_CENTRE_INTERVALS[previous_month]
            end
            source_row = 12 * year + source_month
            target_row = 12 * year + target_month
            fraction = T(offset) / T(denominator)
            @views daily[day, :] .= values[source_row, :] .+
                fraction .* (values[target_row, :] .- values[source_row, :])
            day += 1
        end
    end
    return daily
end

function _daily_deposition_forcing(values::AbstractMatrix{T}, daily_days::Integer) where {T <: AbstractFloat}
    size(values, 1) == daily_days && return values
    return _lpjml_monthly_to_daily(values, daily_days)
end

"""
ClimateDataLoader(climate, data_index, device)

Extract climate forcing slices for selected grid points and years. Optional
NO3/NH4 deposition may be daily or monthly; monthly inputs are expanded
internally with LPJmL's no-leap, month-centred interpolation.
"""
function ClimateDataLoader(climate::NamedTuple, 
                           data_index::Vector{Int},
                           device;
                           T::Type{<:AbstractFloat} = Float32,
)

    loaded_climate = (
        temp_spinup = T.(climate.temp_spinup[:, data_index]),
        temp = T.(climate.temp[:, data_index]),
        prec = T.(climate.prec[:, data_index]),
        sw = T.(climate.swdown[:, data_index]),
        lw = T.(climate.lwnet[:, data_index]),
        co2 = T.(climate.co2),
    )

    # Normalize optional wind forcing to one internal field name while keeping
    # existing climate archives (which do not contain wind) fully compatible.
    if hasproperty(climate, :windspeed)
        loaded_climate = merge(
            loaded_climate,
            (wind = T.(climate.windspeed[:, data_index]),),
        )
    end

    for name in (:no3_deposition, :nh4_deposition)
        hasproperty(climate, name) || continue
        values = T.(getproperty(climate, name)[:, data_index])
        loaded_climate = merge(
            loaded_climate,
            (; name => _daily_deposition_forcing(values, size(loaded_climate.temp, 1))),
        )
    end

    if hasproperty(climate, :co2_daily)
        loaded_climate = merge(loaded_climate, (co2_daily = climate.co2_daily,))
    end

    return _adapt_to_device(device, loaded_climate)
end
