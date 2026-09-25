"""
    dynamic_sowing_date!(crop, climbuf, cft, day; prescribed_winter_type=nothing)

Resolve today's sowing event for `sowing_mode=:dynamic_sdate` without changing
the prescribed-sowing path. An LPJmL-style climate calendar is updated after
each complete climate year from rolling monthly temperature, precipitation, and
potential-evaporation climatologies. The calendar supplies a climate candidate
and its first valid value is retained as a reference. The resulting sowing date
is the prescribed spatial date plus the bounded circular difference between
the current and reference climate candidates. This preserves the local crop
calendar while allowing climate-driven shifts.

The CFT method is also respected: `TEMP_WTYP_CALC_SDATE` is used for wheat,
`TEMP_PREC_CALC_SDATE` for maize, and `PREC_CALC_SDATE` for rice and soybean.
As in LPJmL's current `NO_FIXED_SDATE` path, the seasonality class chooses
the temperature or wet-season trigger; the CFT method only selects the winter
or spring temperature threshold when a temperature trigger is active.

Until a full temperature/precipitation/PET climatology exists, the prescribed
sowing month is used as a deterministic bootstrap month.
"""
function dynamic_sowing_date!(
    state::ModelState,
    climbuf::ClimBuf,
    cft::CFTParameters,
    day::Integer;
    irrigated::Bool = false,
    prescribed_winter_type = nothing,
)
    calendar = crop_calendar_input(state)
    return _dynamic_sowing_date!(
        calendar.sowing_date, calendar.prescribed_sowing_date,
        crop_prognostic(state).phenology.is_growing, climbuf, Int(day),
    )
end

function dynamic_sowing_date!(
    crop::Crop,
    climbuf::ClimBuf,
    cft::CFTParameters,
    day::Integer;
    irrigated::Bool = false,
    prescribed_winter_type = nothing,
)
    calendar = crop.auxiliary.calendar
    return _dynamic_sowing_date!(
        calendar.sowing_date, calendar.prescribed_sowing_date,
        crop.state.phenology.is_growing, climbuf, Int(day),
    )
end

"""
    update_dynamic_sowing_calendar!(climbuf, cft, winter_type)

Update the LPJmL climate calendar once after a complete climate year. The
resulting seasonality class and sowing month are persistent climate state, not
daily event state, so daily sowing cannot reclassify itself mid-season.
"""
function update_dynamic_sowing_calendar!(climbuf::ClimBuf,
                                          cft::CFTParameters,
                                          winter_type::AbstractVector{Bool})
    launch_1D!(
        dynamic_sowing_calendar_kernel!,
        climbuf.seasonality_type,
        climbuf.sowing_month,
        climbuf.climate_sowing_day,
        climbuf.reference_sowing_day,
        climbuf.mtemp20,
        climbuf.mprec20,
        climbuf.mpet20,
        winter_type,
        cft,
    )
    return nothing
end

function _dynamic_sowing_date!(
    sowing_date,
    prescribed_sowing_date,
    is_growing,
    climbuf::ClimBuf,
    day::Int,
)
    launch_1D!(
        dynamic_sowing_date_kernel!,
        sowing_date,
        prescribed_sowing_date,
        is_growing,
        climbuf.climate_sowing_day,
        climbuf.reference_sowing_day,
        day,
    )
    return nothing
end

const _DYNAMIC_NDAYMONTH = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
const _DYNAMIC_MONTH_STARTS = (1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
# LPJmL's no-leap interpolation uses the distance between adjacent monthly
# centres, rather than the number of days in the month itself.
const _DYNAMIC_INTERPOLATION_DENOM = (30, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
const _DYNAMIC_SOWING_MAXIMUM_SHIFT_DAYS = 30

@inline function _dynamic_sowing_month_from_day(day::Integer)
    month = 12
    for candidate in 1:12
        if day >= _DYNAMIC_MONTH_STARTS[candidate]
            month = candidate
        end
    end
    return month
end

@inline function _dynamic_sowing_circular_difference(day::Integer, reference::Integer)
    difference = Int(day) - Int(reference)
    difference > 182 && return difference - 365
    difference < -182 && return difference + 365
    return difference
end

@inline function _dynamic_sowing_shifted_day(
    prescribed::Integer,
    climate_day::Integer,
    reference_day::Integer,
)
    shift = clamp(
        _dynamic_sowing_circular_difference(climate_day, reference_day),
        -_DYNAMIC_SOWING_MAXIMUM_SHIFT_DAYS,
        _DYNAMIC_SOWING_MAXIMUM_SHIFT_DAYS,
    )
    day = Int(prescribed) + shift
    day < 1 && return day + 365
    day > 365 && return day - 365
    return day
end

@inline function _dynamic_sowing_month_from_precipitation(mprec20, mpet20, cell, fallback)
    valid = true
    for month in 1:12
        valid &= mprec20[month, cell] > -9998
        valid &= mpet20[month, cell] > -9998
    end
    valid || return fallback
    best_month = 1
    best_total = zero(eltype(mprec20))
    for month in 1:12
        total = zero(eltype(mprec20))
        for offset in 0:3
            index = mod1(month + offset, 12)
            pet = mpet20[index, cell]
            total += pet > zero(pet) ? mprec20[index, cell] / pet : zero(pet)
        end
        if total > best_total
            best_total = total
            best_month = month
        end
    end
    return best_month
end

@inline function _dynamic_sowing_interpolated_temperature(
    mtemp20::AbstractMatrix{T}, cell, month, day,
) where {T <: AbstractFloat}
    # This mirrors LPJmL's `interpolate(mtemp20, month, day)`: monthly
    # climatologies are centred on the middle of each month, so the value at
    # the end of a month is interpolated from the two adjacent monthly means.
    # The unusual 30/29 denominators are LPJmL's no-leap-calendar spacing
    # between adjacent monthly centres and are kept deliberately here.
    current_month = month
    if day >= _DYNAMIC_NDAYMONTH[month] ÷ 2
        day -= _DYNAMIC_NDAYMONTH[month] ÷ 2
        next_month = month == 12 ? 1 : month + 1
    else
        next_month = month
        current_month = month == 1 ? 12 : month - 1
        day += (_DYNAMIC_NDAYMONTH[current_month] + 1) ÷ 2
    end
    current = mtemp20[current_month, cell]
    next = mtemp20[next_month, cell]
    return current + T(day) * (next - current) /
        T(_DYNAMIC_INTERPOLATION_DENOM[current_month])
end

@inline function _dynamic_sowing_month_from_temperature(mtemp20, cell, threshold, cooling, fallback)
    valid = true
    for month in 1:12
        valid &= mtemp20[month, cell] > -9998
    end
    valid || return fallback
    for month in 1:12
        previous_month = mod1(month - 1, 12)
        previous = mtemp20[previous_month, cell]
        current = mtemp20[month, cell]
        if cooling ? (current < threshold && previous >= threshold) :
                     (current > threshold && previous <= threshold)
            # LPJmL assigns the transition to the previous month when the
            # interpolated value at that month's end has already crossed the
            # threshold; otherwise it assigns the current month.
            previous_end = _dynamic_sowing_interpolated_temperature(
                mtemp20, cell, previous_month,
                _DYNAMIC_NDAYMONTH[previous_month],
            )
            return cooling ?
                (previous_end < threshold ? previous_month : month) :
                (previous_end > threshold ? previous_month : month)
        end
    end
    return 1
end

@inline function _dynamic_sowing_day_from_temperature(
    mtemp20::AbstractMatrix{T},
    cell,
    threshold,
    cooling,
    fallback,
) where {T <: AbstractFloat}
    valid = true
    for month in 1:12
        valid &= mtemp20[month, cell] > T(-9998)
    end
    valid || return fallback

    previous_day = 365
    previous_month = _dynamic_sowing_month_from_day(previous_day)
    previous_temperature = _dynamic_sowing_interpolated_temperature(
        mtemp20, cell, previous_month,
        previous_day - _DYNAMIC_MONTH_STARTS[previous_month] + 1,
    )
    for day in 1:365
        month = _dynamic_sowing_month_from_day(day)
        dayofmonth = day - _DYNAMIC_MONTH_STARTS[month] + 1
        temperature = _dynamic_sowing_interpolated_temperature(
            mtemp20, cell, month, dayofmonth,
        )
        crossed = cooling ?
            temperature < threshold && previous_temperature >= threshold :
            temperature > threshold && previous_temperature <= threshold
        crossed && return day
        previous_temperature = temperature
    end
    return fallback
end

const _DYNAMIC_NO_SEASONALITY = Int32(0)
const _DYNAMIC_PRECIPITATION = Int32(1)
const _DYNAMIC_PRECIPITATION_TEMPERATURE = Int32(2)
const _DYNAMIC_TEMPERATURE = Int32(3)
const _DYNAMIC_TEMPERATURE_PRECIPITATION = Int32(4)
const _DYNAMIC_CLIMATOLOGY_UNAVAILABLE = Int32(-1)

@inline function _dynamic_sowing_seasonality(mtemp20, mprec20, mpet20, cell)
    T = eltype(mtemp20)
    temp_mean = zero(T)
    precip_mean = zero(T)
    temp_valid = true
    precip_valid = true
    pet_valid = true
    for month in 1:12
        temperature = mtemp20[month, cell]
        temp_valid &= temperature > T(-9998)
        temp_mean += temperature + T(273.15)
        precipitation = mprec20[month, cell]
        precip_valid &= precipitation > T(-9998)
        precip_mean += precipitation
        pet_valid &= mpet20[month, cell] > T(-9998)
    end
    # Keep the bootstrap state distinct from a valid, non-seasonal climate.
    # LPJmL's initialized buffers are unavailable until their first complete
    # climatology update; the deterministic prescribed-month fallback is only
    # for that bootstrap period.
    temp_valid && precip_valid && pet_valid || return _DYNAMIC_CLIMATOLOGY_UNAVAILABLE
    temp_mean /= T(12)
    precip_mean /= T(12)
    temp_mean > zero(T) && precip_mean > zero(T) || return _DYNAMIC_NO_SEASONALITY

    temp_variance = zero(T)
    precip_variance = zero(T)
    coldest = typemax(T)
    for month in 1:12
        temperature = mtemp20[month, cell] + T(273.15)
        precipitation = mprec20[month, cell]
        coldest = min(coldest, temperature - T(273.15))
        temp_variance += (temperature - temp_mean)^2
        precip_variance += (precipitation - precip_mean)^2
    end
    # LPJmL uses the sample coefficient of variation and the fixed thresholds
    # var_temp > 0.010, var_prec > 0.4, coldest month <= 10 °C.
    cv_temp = sqrt(temp_variance / T(11)) / temp_mean
    cv_precip = sqrt(precip_variance / T(11)) / precip_mean
    temp_seasonal = cv_temp > T(0.010)
    precip_seasonal = cv_precip > T(0.4)
    if precip_seasonal && !temp_seasonal
        return _DYNAMIC_PRECIPITATION
    elseif precip_seasonal && temp_seasonal && coldest > T(10)
        return _DYNAMIC_PRECIPITATION_TEMPERATURE
    elseif temp_seasonal && !precip_seasonal
        return _DYNAMIC_TEMPERATURE
    elseif temp_seasonal && precip_seasonal
        return _DYNAMIC_TEMPERATURE_PRECIPITATION
    end
    return _DYNAMIC_NO_SEASONALITY
end

@kernel inbounds = true function dynamic_sowing_calendar_kernel!(
    seasonality_type::AbstractVector{Int32},
    sowing_month::AbstractVector{Int32},
    climate_sowing_day::AbstractVector{Int32},
    reference_sowing_day::AbstractVector{Int32},
    mtemp20::AbstractMatrix{T},
    mprec20::AbstractMatrix{T},
    mpet20::AbstractMatrix{T},
    winter_type::AbstractVector{Bool},
    cft::CFTParameters,
) where {T <: AbstractFloat}
    cell = @index(Global)
    seasonality = _dynamic_sowing_seasonality(mtemp20, mprec20, mpet20, cell)
    seasonality_type[cell] = seasonality
    climate_day = 0
    if seasonality == _DYNAMIC_CLIMATOLOGY_UNAVAILABLE
        sowing_month[cell] = Int32(0)
    else
        method = cft.sowing_date.method
        if method == SDATE_NO_CALC || method == SDATE_MULTICROP
            sowing_month[cell] = Int32(0)
        elseif seasonality == _DYNAMIC_NO_SEASONALITY
            sowing_month[cell] = Int32(1)
            climate_day = 1
        elseif seasonality == _DYNAMIC_PRECIPITATION ||
               seasonality == _DYNAMIC_PRECIPITATION_TEMPERATURE
            sowing_month[cell] = Int32(_dynamic_sowing_month_from_precipitation(
                mprec20, mpet20, cell, 1,
            ))
            climate_day = _DYNAMIC_MONTH_STARTS[sowing_month[cell]]
        else
            winter_sowing = method == SDATE_TEMPERATURE_WINTER && winter_type[cell]
            threshold = winter_sowing ? cft.sowing_date.temp_fall : cft.sowing_date.temp_spring
            sowing_month[cell] = Int32(_dynamic_sowing_month_from_temperature(
                mtemp20, cell, threshold, winter_sowing, 1,
            ))
            climate_day = _dynamic_sowing_day_from_temperature(
                mtemp20, cell, threshold, winter_sowing,
                _DYNAMIC_MONTH_STARTS[sowing_month[cell]],
            )
        end
    end
    climate_sowing_day[cell] = Int32(climate_day)
    if climate_day > 0 && reference_sowing_day[cell] == 0
        reference_sowing_day[cell] = Int32(climate_day)
    end
end

@kernel inbounds = true function dynamic_sowing_date_kernel!(
    sowing_date::AbstractVector{I},
    prescribed_sowing_date::AbstractVector{I},
    is_growing::AbstractVector{G},
    climate_sowing_day::AbstractVector{Int32},
    reference_sowing_day::AbstractVector{Int32},
    day::Int,
) where {I <: Integer, G <: Integer}
    cell = @index(Global)

    # The trigger is ephemeral: cultivate! consumes it today, then the next
    # dynamic resolution clears it. Before a complete climate calendar exists,
    # or for an unsupported calendar method, use the prescribed date exactly.
    sowing_date[cell] = zero(I)
    if iszero(is_growing[cell])
        prescribed = prescribed_sowing_date[cell]
        climate_day = climate_sowing_day[cell]
        reference_day = reference_sowing_day[cell]
        target_day = climate_day > 0 && reference_day > 0 ?
            _dynamic_sowing_shifted_day(prescribed, climate_day, reference_day) :
            Int(prescribed)
        sowing_date[cell] = day == target_day ? I(day) : zero(I)
    end
end
