function spin_up_climbuf!(PFT::PftParameters, 
                          temp_spinup::AbstractArray{T}, 
                          climbuf::ClimBuf,
                          year_spinup::Integer,
                          device
) where {T <: AbstractFloat}
    for i = 1 : year_spinup
        year_temp = temp_spinup[365*(i-1)+1 : 365*i, :]
        for day in axes(year_temp, 1)
            daily_climbuf!(year_temp[day, :], climbuf.temp)
        end
        climbuf.V_req_a .= zero(T)
        annual_climbuf!(year_temp, climbuf, PFT, device)
    end

end

function update_climbuf!(PFT::PftParameters, 
                         temp::AbstractArray{T},
                         climate_temp::AbstractArray{M},
                         climbuf::ClimBuf,
                         day::Integer,
                         device
)where {T <: AbstractFloat, M <: AbstractFloat, }

    daily_climbuf!(temp, climbuf.temp)

    if day > 1 && day % 365 == 1
        year_temp = climate_temp[day-365:day-1, :]
        climbuf.V_req_a .= 0.0f0
        annual_climbuf!(year_temp, climbuf, PFT, device)
    end

end