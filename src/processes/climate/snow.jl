"""
snow!(soil, dailyWeather)

Update snowpack, snow height, and snow cover fraction from daily temperature
and precipitation forcing.
"""
function snow!(soil,
               dailyWeather::DailyWeather;
               snowparams::SnowParams = snowparams,
               lpjmlparams::LPJmLParams = lpjmlparams
)

    kernel_params = kernel_constant(
        dailyWeather.temp, (; snowparams, lpjmlparams),
    )

    launch_1D!(
        snow_kernel!,
        dailyWeather.temp,
        dailyWeather.prec,
        soil_snow_prognostic(soil).pack,
        soil_snow_fluxes(soil).melt,
        soil_snow_fluxes(soil).sublimation,
        soil_snow_fluxes(soil).runoff,
        soil_snow_prognostic(soil).height,
        soil_snow_prognostic(soil).fraction,
        kernel_params
    )

end

function snow!(state::ModelState;
               snowparams::SnowParams = snowparams,
               lpjmlparams::LPJmLParams = lpjmlparams)
    weather = state.inputs.weather
    snow_state = state.prognostic.soil.snow
    snow_fluxes = state.fluxes.soil.snow
    launch_1D!(
        snow_kernel!, weather.temp, weather.prec,
        snow_state.pack, snow_fluxes.melt, snow_fluxes.sublimation,
        snow_fluxes.runoff, snow_state.height, snow_state.fraction,
        kernel_constant(weather.temp, (; snowparams, lpjmlparams)),
    )
    return nothing
end

"""
    compute_snowfall(pack, precipitation, temperature, threshold, maximum_pack)

Partition precipitation at the snow threshold. Returns the updated snowpack,
remaining liquid precipitation, and overflow runoff in millimetres. This is a
cell-local LPJmL rule; it deliberately has no array access or state mutation.
"""
@inline function compute_snowfall(pack::T,
                                  precipitation::T,
                                  temperature::T,
                                  threshold::T,
                                  maximum_pack::T) where {T <: AbstractFloat}
    temperature < threshold || return pack, precipitation, zero(T)
    updated_pack = pack + precipitation
    runoff = max(updated_pack - maximum_pack, zero(T))
    return min(updated_pack, maximum_pack), zero(T), runoff
end

"""
    compute_snow_sublimation(pack)

Remove LPJmL's fixed daily 0.1 mm sublimation only when enough snow is
present. Returning both quantities makes the stock--flux relationship explicit.
"""
@inline function compute_snow_sublimation(pack::T) where {T <: AbstractFloat}
    pack > T(0.1) || return pack, zero(T)
    return pack - T(0.1), T(0.1)
end

"""
    compute_snow_melt(pack, temperature, threshold, params)

Compute the conductive snow-skin melt for one daily step. The returned melt is
not clipped to the snowpack because the original LPJmL-compatible update first
applies it and then clears a near-empty pack using `snow_epsilon`.
"""
@inline function compute_snow_melt(pack::T,
                                   temperature::T,
                                   threshold::T,
                                   snow_skin_depth::T,
                                   thermal_diffusivity::T,
                                   conductivity::T,
                                   fusion_capacity::T) where {T <: AbstractFloat}
    (pack > T(1e-7) && temperature > zero(T)) || return zero(T)
    timestep_seconds = T(24 * 3600)
    depth = min(pack, snow_skin_depth)
    temperature_increment = thermal_diffusivity * timestep_seconds / (depth * depth) *
        T(1_000_000) * (temperature - threshold)
    heat_flux = conductivity * (threshold + temperature_increment) / depth * T(1000)
    melt_heat = min(heat_flux * timestep_seconds, depth * T(1e-3) * fusion_capacity)
    return melt_heat / fusion_capacity * T(1000)
end

"""
    compute_snow_geometry(pack, water_to_snow, roughness)

Return the normalized snowpack, physical snow height (m), and fractional ground
cover. The normalization prevents small negative residuals from propagating.
"""
@inline function compute_snow_geometry(pack::T,
                                       water_to_snow::T,
                                       roughness::T) where {T <: AbstractFloat}
    pack > T(1e-7) || return zero(T), zero(T), zero(T)
    height = water_to_snow * (pack / T(1000))
    fraction = height / (height + T(0.5) * roughness)
    return pack, height, fraction
end



@kernel inbounds = true function snow_kernel!(
                              temp::AbstractArray{T},
                              prec::AbstractArray{T},
                              soil_snowpack::AbstractArray{T},
                              soil_snowmelt::AbstractArray{T},
                              soil_snow_sublimation::AbstractArray{T},
                              soil_snow_runoff::AbstractArray{T},
                              soil_snowheight::AbstractArray{T},
                              soil_snowfraction::AbstractArray{T},
                              fixed_parameters
) where {T <: AbstractFloat}

    cell = @index(Global)
    kernel_params = kernel_value(fixed_parameters)

    @unpack tsnow, snow_skin_depth, th_diff_snow, lambda_snow, c_water2ice, c_watertosnow, c_roughness= kernel_params.snowparams
    @unpack maxsnowpack = kernel_params.lpjmlparams

    soil_snowmelt[cell] = zero(T)
    soil_snow_sublimation[cell] = zero(T)
    soil_snow_runoff[cell] = zero(T)

    # Phase partition and sublimation are independent cell-local stock updates.
    pack, liquid_precipitation, runoff = compute_snowfall(
        soil_snowpack[cell], prec[cell], temp[cell], T(tsnow), T(maxsnowpack),
    )
    pack, sublimation = compute_snow_sublimation(pack)
    soil_snowpack[cell] = pack
    prec[cell] = liquid_precipitation
    soil_snow_runoff[cell] = runoff
    soil_snow_sublimation[cell] = sublimation

    # snow layer is insulating
    melt = compute_snow_melt(
        soil_snowpack[cell], temp[cell], T(tsnow), T(snow_skin_depth),
        T(th_diff_snow), T(lambda_snow), T(c_water2ice),
    )
    soil_snowmelt[cell] = melt
    soil_snowpack[cell] -= melt

    # Add melt water to rainfall before interception and infiltration.
    prec[cell] += soil_snowmelt[cell]

    # calculate snow height and fraction of snow coverage
    pack, height, fraction = compute_snow_geometry(
        soil_snowpack[cell], T(c_watertosnow), T(c_roughness),
    )
    soil_snowpack[cell] = pack
    soil_snowheight[cell] = height
    soil_snowfraction[cell] = fraction

end
