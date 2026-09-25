"""
spin_up_climbuf!(CFT, climate, climbuf, day, lat, temp, lwnet, swdown, output)

Run climate-buffer spin-up for one step before full crop process integration.
"""
function spin_up_climbuf!(CFT::CFTParameters,
                          temp_spinup::AbstractArray{T}, 
                          climbuf::ClimBuf;
                          prec_spinup = nothing,
                          year_spinup = 1
) where {T <: AbstractFloat}
    prec_spinup === nothing || size(prec_spinup) == size(temp_spinup) || throw(DimensionMismatch(
        "precipitation spin-up must match temperature spin-up dimensions",
    ))
    for i = 1 : year_spinup
        year_temp = temp_spinup[365*(i-1)+1 : 365*i, :]
        for day in axes(year_temp, 1)
            daily_climbuf!(year_temp[day, :], climbuf.temp)
        end
        year_prec = prec_spinup === nothing ? nothing : prec_spinup[365*(i-1)+1 : 365*i, :]
        annual_climbuf!(year_temp, climbuf, CFT; daily_prec = year_prec)
    end

end


"""
update_climbuf!(CFT, climbuf, day, lat, temp, lwnet, swdown)

Update climate-buffer and PET diagnostics during daily simulation.
"""
function update_climbuf!(CFT::CFTParameters,
                         temp::AbstractArray{T},
                         climbuf::ClimBuf,
                         day::Integer;
                         prec = nothing,
                         dynamic_sowing::Bool = false,
                         winter_type = nothing,
                         update_vernalization_requirement::Bool = true,
)where {T <: AbstractFloat}

    if day > 1 && day % 365 == 1
        annual_climbuf!(
            climbuf.atemp, climbuf, CFT;
            daily_prec = climbuf.prec,
            daily_pet = climbuf.apet,
            update_vernalization_requirement,
        )
        if dynamic_sowing
            isnothing(winter_type) && throw(ArgumentError(
                "dynamic sowing requires winter-type input",
            ))
            update_dynamic_sowing_calendar!(climbuf, CFT, winter_type)
        end
    end
    
    if day % 365 == 0
        day = 365
    else
        day = day % 365
    end
    
    daily_climbuf!(
        temp, climbuf.temp;
        annual_temperature = climbuf.atemp,
        annual_day = day,
    )
    prec === nothing || daily_climbuf!(
        prec, climbuf.prec;
        annual_temperature = climbuf.prec,
        annual_day = day,
    )

end

"""
    record_potential_evaporation!(climbuf, equilibrium_evaporation, day, params)

Record LPJmL-style potential evaporation (`eeq * ALPHAM`) for the monthly
dynamic-sowing climatology. This diagnostic does not alter radiation,
transpiration, or soil-water calculations.
"""
function record_potential_evaporation!(climbuf::ClimBuf,
                                       equilibrium_evaporation::AbstractArray{T},
                                       day::Integer,
                                       params::LPJmLParams = lpjmlparams,
) where {T <: AbstractFloat}
    annual_day = day % 365 == 0 ? 365 : day % 365
    daily_climbuf!(
        equilibrium_evaporation, climbuf.apet;
        annual_temperature = climbuf.apet,
        annual_day,
        shift_history = false,
        scale = T(params.ALPHAM),
    )
    return nothing
end
