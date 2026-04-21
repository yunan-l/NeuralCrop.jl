function readclimate!(climate::NamedTuple,
                      dailyWeather::DailyWeather,
                      day::Integer
)

    dailyWeather.temp = climate.temp[day, :]
    dailyWeather.prec = climate.prec[day, :]
    dailyWeather.swr = climate.sw[day, :]
    dailyWeather.lwr = climate.lw[day, :]
    dailyWeather.temp_n = climate.temp_n[day, :]
    # dailyWeather.prec_n = climate.prec_n[day, :]
    dailyWeather.swr_n = climate.sw_n[day, :]
    dailyWeather.lwr_n = climate.lw_n[day, :]
    dailyWeather.annual_co2 = ppm2Pa(climate.co2[[div(day-1, 365) + 1]])

end

function readclimate!(climate::NamedTuple,
                      dailyWeather::DailyWeather,
                      CO2::AbstractArray{T}, 
                      day::Integer
) where {T <: AbstractFloat}

    dailyWeather.temp = climate.temp[day, :]
    dailyWeather.prec = climate.prec[day, :]
    dailyWeather.swr = climate.sw[day, :]
    dailyWeather.lwr = climate.lw[day, :]
    dailyWeather.temp_n = climate.temp_n[day, :]
    # dailyWeather.prec_n = climate.prec_n[day, :]
    dailyWeather.swr_n = climate.sw_n[day, :]
    dailyWeather.lwr_n = climate.lw_n[day, :]
    dailyWeather.daily_co2 = ppm2Pa(CO2[day, :])
    
end