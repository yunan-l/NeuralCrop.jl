"""
    readclimate!(climate, dailyWeather, day)

Read one day of climate forcing and return the active CO₂ buffer. A vector
`climate.co2` is interpreted as a global annual series shared by all cells; a
matrix is interpreted as daily forcing with shape `(day, cell)`.
"""

function readclimate!(climate::NamedTuple,
                      dailyWeather::DailyWeather,
                      day::Integer)
    has_wind = hasproperty(climate, :wind)
    has_no3_deposition = hasproperty(climate, :no3_deposition)
    has_nh4_deposition = hasproperty(climate, :nh4_deposition)
    co2_daily = hasproperty(climate, :co2_daily) && climate.co2_daily
    wind = has_wind ? climate.wind : climate.temp
    no3_deposition = has_no3_deposition ? climate.no3_deposition : climate.temp
    nh4_deposition = has_nh4_deposition ? climate.nh4_deposition : climate.temp
    default_wind = eltype(dailyWeather.temp)(lpjmlparams.volatil_wind)
    if ndims(climate.co2) == 1
        launch_1D!(
            read_annual_climate_kernel!,
            dailyWeather.temp,
            dailyWeather.prec,
            dailyWeather.swr,
            dailyWeather.lwr,
            dailyWeather.wind,
            dailyWeather.no3_deposition,
            dailyWeather.nh4_deposition,
            dailyWeather.annual_co2,
            climate.temp,
            climate.prec,
            climate.sw,
            climate.lw,
            wind,
            no3_deposition,
            nh4_deposition,
            climate.co2,
            kernel_constant(
                dailyWeather.temp,
                (;
                    day,
                    has_wind,
                    has_no3_deposition,
                    has_nh4_deposition,
                    co2_daily,
                    default_wind,
                ),
            ),
        )
        return dailyWeather.annual_co2
    elseif ndims(climate.co2) == 2
        launch_1D!(
            read_daily_climate_kernel!,
            dailyWeather.temp,
            dailyWeather.prec,
            dailyWeather.swr,
            dailyWeather.lwr,
            dailyWeather.wind,
            dailyWeather.no3_deposition,
            dailyWeather.nh4_deposition,
            dailyWeather.daily_co2,
            climate.temp,
            climate.prec,
            climate.sw,
            climate.lw,
            wind,
            no3_deposition,
            nh4_deposition,
            climate.co2,
            kernel_constant(
                dailyWeather.temp,
                (; day, has_wind, has_no3_deposition, has_nh4_deposition, default_wind),
            ),
        )
        return dailyWeather.daily_co2
    else
        throw(ArgumentError("climate.co2 must be a vector or a (day, cell) matrix"))
    end
end

@kernel inbounds = true function read_annual_climate_kernel!(
    temperature::AbstractVector{T},
    precipitation::AbstractVector{T},
    shortwave::AbstractVector{T},
    longwave::AbstractVector{T},
    wind::AbstractVector{T},
    no3_deposition::AbstractVector{T},
    nh4_deposition::AbstractVector{T},
    annual_co2::AbstractVector{T},
    temperature_forcing::AbstractMatrix{T},
    precipitation_forcing::AbstractMatrix{T},
    shortwave_forcing::AbstractMatrix{T},
    longwave_forcing::AbstractMatrix{T},
    wind_forcing::AbstractMatrix{T},
    no3_deposition_forcing::AbstractMatrix{T},
    nh4_deposition_forcing::AbstractMatrix{T},
    co2_forcing::AbstractVector{T},
    fixed_parameters,
) where {T <: AbstractFloat}
    cell = @index(Global)
    parameters = kernel_value(fixed_parameters)
    day = parameters.day
    has_wind = parameters.has_wind
    has_no3_deposition = parameters.has_no3_deposition
    has_nh4_deposition = parameters.has_nh4_deposition
    co2_daily = parameters.co2_daily
    default_wind = T(parameters.default_wind)
    temperature[cell] = temperature_forcing[day, cell]
    precipitation[cell] = precipitation_forcing[day, cell]
    shortwave[cell] = shortwave_forcing[day, cell]
    longwave[cell] = longwave_forcing[day, cell]
    wind[cell] = has_wind ? wind_forcing[day, cell] : default_wind
    no3_deposition[cell] = has_no3_deposition ? no3_deposition_forcing[day, cell] : zero(T)
    nh4_deposition[cell] = has_nh4_deposition ? nh4_deposition_forcing[day, cell] : zero(T)
    if cell == 1
        co2_index = co2_daily ? day : div(day - 1, 365) + 1
        annual_co2[1] = co2_forcing[co2_index] * T(0.1)
    end
end

@kernel inbounds = true function read_daily_climate_kernel!(
    temperature::AbstractVector{T},
    precipitation::AbstractVector{T},
    shortwave::AbstractVector{T},
    longwave::AbstractVector{T},
    wind::AbstractVector{T},
    no3_deposition::AbstractVector{T},
    nh4_deposition::AbstractVector{T},
    daily_co2::AbstractVector{T},
    temperature_forcing::AbstractMatrix{T},
    precipitation_forcing::AbstractMatrix{T},
    shortwave_forcing::AbstractMatrix{T},
    longwave_forcing::AbstractMatrix{T},
    wind_forcing::AbstractMatrix{T},
    no3_deposition_forcing::AbstractMatrix{T},
    nh4_deposition_forcing::AbstractMatrix{T},
    co2_forcing::AbstractMatrix{T},
    fixed_parameters,
) where {T <: AbstractFloat}
    cell = @index(Global)
    parameters = kernel_value(fixed_parameters)
    day = parameters.day
    has_wind = parameters.has_wind
    has_no3_deposition = parameters.has_no3_deposition
    has_nh4_deposition = parameters.has_nh4_deposition
    default_wind = T(parameters.default_wind)
    temperature[cell] = temperature_forcing[day, cell]
    precipitation[cell] = precipitation_forcing[day, cell]
    shortwave[cell] = shortwave_forcing[day, cell]
    longwave[cell] = longwave_forcing[day, cell]
    wind[cell] = has_wind ? wind_forcing[day, cell] : default_wind
    no3_deposition[cell] = has_no3_deposition ? no3_deposition_forcing[day, cell] : zero(T)
    nh4_deposition[cell] = has_nh4_deposition ? nh4_deposition_forcing[day, cell] : zero(T)
    daily_co2[cell] = co2_forcing[day, cell] * T(0.1)
end
