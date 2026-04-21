function snow!(soil::Soil,
               dailyWeather::DailyWeather;
               snowparams::SnowParams = snowparams,
               lpjmlparams::LPJmLParams = lpjmlparams,
)
    backend = KernelAbstractions.get_backend(dailyWeather.temp)
    
    kernel = snow_kernel!(backend)
    
    kernel(snowparams, lpjmlparams, dailyWeather.temp, dailyWeather.prec, soil.snowpack, soil.snowheight, soil.snowfraction, ndrange=length(dailyWeather.temp))
    
    KernelAbstractions.synchronize(backend)

end



@kernel function snow_kernel!(snowparams::SnowParams,
                              lpjmlparams::LPJmLParams,
                              temp::AbstractArray{T},
                              prec::AbstractArray{T},
                              soil_snowpack::AbstractArray{T},
                              soil_snowheight::AbstractArray{T},
                              soil_snowfraction::AbstractArray{T},
) where {T <: AbstractFloat}
    
    cell = @index(Global)

    @unpack tsnow, snow_skin_depth, th_diff_snow, lambda_snow, c_water2ice, c_watertosnow, c_roughness= snowparams
    @unpack maxsnowpack = lpjmlparams

    # precipitation falls as snow
    if temp[cell] < tsnow
        soil_snowpack[cell] += prec[cell]
        if soil_snowpack[cell] > maxsnowpack
            # runoff[cell] = soil_snowpack[cell] - maxsnowpack
            soil_snowpack[cell] = maxsnowpack
        end
        prec[cell] = zero(T)
    end

    # sublimation of snow
    if soil_snowpack[cell] > T(0.1)
        soil_snowpack[cell] -= T(0.1)
        # evap[cell] += T(0.1)
    # else
        # evap[cell] = T(0.1)
    
    end

    # snow layer is insulating
    timestep2sec = T(24.0 * 3600.0)
    if soil_snowpack[cell] > T(1e-7)
        if temp[cell] > zero(T)
            depth = min(soil_snowpack[cell], snow_skin_depth)
            dT = th_diff_snow * timestep2sec / (depth * depth) * T(1000000.0) * (temp[cell] - tsnow)
            heatflux = lambda_snow * (tsnow - zero(T) + dT) / depth * T(1000)
            melt_heat = min(heatflux * timestep2sec, depth * T(1e-3) * c_water2ice) #[J/m2]
            # snowmelt[cell] += melt_heat / c_water2ice * T(1000)
            soil_snowpack[cell] -= melt_heat / c_water2ice * T(1000)
            if soil_snowpack[cell] < T(1e-7)
                soil_snowpack[cell] = zero(T)
                soil_snowheight[cell] = zero(T)
                soil_snowfraction[cell] = zero(T)
            end
        end
    end

    # calculate snow height and fraction of snow coverage 
    if soil_snowpack[cell] > T(1e-7)
        HS = c_watertosnow * (soil_snowpack[cell] / T(1000.0)) # mm -> m */
        frsg = HS / (HS+ T(0.5) * c_roughness)
        soil_snowheight[cell] = HS
        soil_snowfraction[cell] = frsg
    else
        soil_snowpack[cell] = zero(T)
        soil_snowheight[cell] = zero(T)
        soil_snowfraction[cell] = zero(T)
    end

end