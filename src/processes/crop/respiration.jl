"""
respiration!(crop, CFT, temp, assim; lpjmlparams=lpjmlparams)

Compute maintenance and growth respiration and update `crop_fluxes(crop).carbon.respiration`.
"""

"""Allocation-free daily respiration using one cell-local CPU/GPU kernel."""
function respiration!(crop,
                      CFT::CFTParameters,
                      air_temperature::AbstractVector{T},
                      soil_temperature::AbstractMatrix{T},
                      gross_assimilation::AbstractArray{T},
                      leaf_respiration::AbstractArray{T};
                      crop_resp_fix::Bool = false,
                      lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat}
    launch_1D!(
        respiration_kernel!,
        crop_fluxes(crop).carbon.respiration,
        crop_prognostic(crop).carbon.root,
        crop_prognostic(crop).carbon.storage,
        crop_prognostic(crop).carbon.pool,
        crop_prognostic(crop).nitrogen.root,
        crop_prognostic(crop).nitrogen.storage,
        crop_prognostic(crop).nitrogen.pool,
        crop_prognostic(crop).phenology.is_growing,
        air_temperature,
        soil_temperature,
        gross_assimilation,
        leaf_respiration,
        kernel_constant(
            crop_fluxes(crop).carbon.respiration, (; CFT, lpjmlparams),
        ),
        crop_resp_fix,
    )
    return nothing
end

# Backward-compatible entry for callers without an explicit soil-temperature
# profile. Daily simulations use the method above.
function respiration!(crop,
                      CFT::CFTParameters,
                      air_temperature::AbstractVector{T},
                      gross_assimilation::AbstractArray{T},
                      leaf_respiration::AbstractArray{T};
                      crop_resp_fix::Bool = false,
                      lpjmlparams::LPJmLParams = lpjmlparams) where {T <: AbstractFloat}
    return respiration!(
        crop, CFT, air_temperature, reshape(air_temperature, 1, :),
        gross_assimilation, leaf_respiration; crop_resp_fix, lpjmlparams = lpjmlparams,
    )
end

"""
    compute_respiration_temperature_response(temperature, e0, response_temperature)

LPJmL's capped Lloyd--Taylor temperature multiplier for maintenance
respiration. Temperatures below -15 °C contribute zero; values above 40 °C
retain the 40 °C response.
"""
@inline function compute_respiration_temperature_response(temperature::T,
                                                           e0::T,
                                                           response_temperature::T) where {T <: AbstractFloat}
    temperature >= T(-15) || return zero(T)
    bounded_temperature = min(temperature, T(40))
    return exp(e0 * (one(T) / (response_temperature + T(10)) -
                     one(T) / (bounded_temperature + response_temperature)))
end

"""
    compute_maintenance_respiration(carbon, coefficient, rate, nitrogen_ratio,
                                    temperature_response)

Compute one organ's daily maintenance respiration. The calling kernel keeps
the separate root/air temperature choices explicit.
"""
@inline function compute_maintenance_respiration(carbon::T,
                                                  coefficient::T,
                                                  rate::T,
                                                  nitrogen_ratio::T,
                                                  temperature_response::T) where {T <: AbstractFloat}
    return carbon * coefficient * rate * nitrogen_ratio * temperature_response
end

"""Return LPJmL-style organ N:C for crop maintenance respiration."""
@inline function compute_respiration_nitrogen_ratio(
    carbon::T,
    nitrogen::T,
    fixed_ratio::T,
    maximum_ratio::T,
    crop_resp_fix::Bool,
) where {T <: AbstractFloat}
    ratio = !crop_resp_fix && carbon > eps(T) ? nitrogen / carbon : fixed_ratio
    return min(max(zero(T), ratio), maximum_ratio)
end

"""
    compute_growth_respiration(assimilation, root, storage, pool, growth_fraction)

Allocate growth respiration only from the positive carbon remainder after leaf
dark respiration and all maintenance terms have been deducted.
"""
@inline function compute_growth_respiration(assimilation::T,
                                             root::T,
                                             storage::T,
                                             pool::T,
                                             growth_fraction::T) where {T <: AbstractFloat}
    return max(zero(T), (assimilation - root - storage - pool) * growth_fraction)
end

@kernel inbounds = true function respiration_kernel!(
    respiration::AbstractVector{T},
    root_carbon::AbstractVector{T},
    storage_carbon::AbstractVector{T},
    pool_carbon::AbstractVector{T},
    root_nitrogen::AbstractVector{T},
    storage_nitrogen::AbstractVector{T},
    pool_nitrogen::AbstractVector{T},
    is_growing::AbstractVector{I},
    air_temperature::AbstractVector{T},
    soil_temperature::AbstractMatrix{T},
    gross_assimilation::AbstractVector{T},
    leaf_respiration::AbstractVector{T},
    fixed_parameters,
    crop_resp_fix::Bool,
) where {T <: AbstractFloat, I <: Integer}
    cell = @index(Global)
    parameters = kernel_value(fixed_parameters)
    CFT = parameters.CFT
    lpjmlparams = parameters.lpjmlparams
    @unpack respcoeff, nc_ratio, ncleaf, ratio = CFT
    @unpack k, r_growth, e0, temp_response = lpjmlparams

    gtemp_air = compute_respiration_temperature_response(
        air_temperature[cell], T(e0), T(temp_response),
    )
    gtemp_soil = compute_respiration_temperature_response(
        soil_temperature[1, cell], T(e0), T(temp_response),
    )
    root_nc = compute_respiration_nitrogen_ratio(
        root_carbon[cell], root_nitrogen[cell], T(nc_ratio.root),
        T(ncleaf.high) / T(ratio.root), crop_resp_fix,
    )
    storage_nc = compute_respiration_nitrogen_ratio(
        storage_carbon[cell], storage_nitrogen[cell], T(nc_ratio.sto),
        T(ncleaf.high) / T(ratio.sto), crop_resp_fix,
    )
    pool_nc = compute_respiration_nitrogen_ratio(
        pool_carbon[cell], pool_nitrogen[cell], T(nc_ratio.pool),
        T(ncleaf.high) / T(ratio.pool), crop_resp_fix,
    )
    root_respiration = compute_maintenance_respiration(
        root_carbon[cell], T(respcoeff), T(k), root_nc, gtemp_soil,
    )
    storage_respiration = compute_maintenance_respiration(
        storage_carbon[cell], T(respcoeff), T(k), storage_nc, gtemp_air,
    )
    pool_respiration = compute_maintenance_respiration(
        pool_carbon[cell], T(respcoeff), T(k), pool_nc, gtemp_air,
    )
    assimilation = gross_assimilation[cell] - leaf_respiration[cell]
    growth_respiration = compute_growth_respiration(
        assimilation, root_respiration, storage_respiration, pool_respiration,
        T(r_growth),
    )
    active = T(is_growing[cell])
    respiration[cell] =
        (root_respiration + storage_respiration + pool_respiration + growth_respiration) * active
end
