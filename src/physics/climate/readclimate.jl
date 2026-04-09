function readclimate!(climate::NamedTuple,
                      day::Integer
)

    temp = climate.temp[day, :]
    prec = climate.prec[day, :]
    swr = climate.sw[day, :]
    lwr = climate.lw[day, :]
    temp_n = climate.temp_n[day, :]
    # prec_n = climate.prec_n[day, :]
    swr_n = climate.sw_n[day, :]
    lwr_n = climate.lw_n[day, :]
    co2 = ppm2Pa(climate.co2[[div(day-1, 365) + 1]])

    return temp, prec, swr, lwr, temp_n, swr_n, lwr_n, co2

end

function readclimate!(climate::NamedTuple,
                      co2::AbstractArray{T}, 
                      day::Integer
) where {T <: AbstractFloat}

    temp = climate.temp[day, :]
    prec = climate.prec[day, :]
    swr = climate.sw[day, :]
    lwr = climate.lw[day, :]
    temp_n = climate.temp_n[day, :]
    prec_n = climate.prec_n[day, :]
    swr_n = climate.sw_n[day, :]
    lwr_n = climate.lw_n[day, :]
    co2 = ppm2Pa(co2[day, :])

    return temp, prec, swr, lwr, temp_n, prec_n, swr_n, lwr_n, co2

end