# using CUDA
"""
photosynthesis_C3!(CFT, photos, crop, pet, co2, temp)

Compute C3 photosynthesis rates and related diagnostic variables.
"""


"""
photosynthesis_C4!(CFT, photos, crop, pet, co2, temp)

Compute C4 photosynthesis rates and related diagnostic variables.
"""

"""Cell-local C3 photosynthesis kernel with no intermediate device arrays."""
function photosynthesis_C3!(CFT::CFTParameters,
                            crop,
                            apar::AbstractArray{T},
                            pet_daylength::AbstractArray{T},
                            temp::AbstractArray{T},
                            co2::AbstractArray{T};
                            lpjmlparams::LPJmLParams = lpjmlparams,
                            photoparams::PhotoParams = photoparams,
                            comp_vcmax = false
) where {T <: AbstractFloat}
    launch_1D!(
        photosynthesis_c3_kernel!,
        crop_fluxes(crop).carbon.gross_assimilation,
        crop_fluxes(crop).carbon.net_assimilation,
        crop_fluxes(crop).carbon.water_limited_assimilation,
        crop_fluxes(crop).carbon.leaf_respiration,
        crop_photosynthesis_auxiliary(crop).potential_vcmax,
        crop_photosynthesis_auxiliary(crop).vcmax,
        crop_photosynthesis_auxiliary(crop).nitrogen_limitation,
        crop_photosynthesis_auxiliary(crop).lambda,
        crop_photosynthesis_auxiliary(crop).temperature_stress,
        apar,
        pet_daylength,
        temp,
        co2,
        kernel_constant(
            crop_fluxes(crop).carbon.gross_assimilation,
            (; CFT, lpjmlparams, photoparams),
        ),
        comp_vcmax,
    )
    return nothing
end

@inline function compute_co_limited_assimilation(
    light_limited::T,
    rubisco_limited::T,
    curvature::T,
    daylength::T,
) where {T <: AbstractFloat}
    discriminant = max(
        zero(T),
        (light_limited + rubisco_limited) * (light_limited + rubisco_limited) -
        T(4) * curvature * light_limited * rubisco_limited,
    )
    return (light_limited + rubisco_limited - sqrt(discriminant)) /
           (T(2) * curvature) * daylength
end

@inline function compute_net_assimilation(
    gross::T,
    leaf_respiration::T,
    daylength::T,
) where {T <: AbstractFloat}
    daily_net = gross - hour2day(daylength) * leaf_respiration
    return max(zero(T), daily_net), daily_net
end

@inline function compute_water_limited_assimilation(
    daily_net::T,
    carbon_mass::T,
    temperature::T,
    pressure::T,
) where {T <: AbstractFloat}
    daily_net <= zero(T) && return zero(T)
    return daily_net / carbon_mass * T(8.314) * (temperature + T(273.15)) /
           pressure * T(1000)
end

@kernel inbounds = true function photosynthesis_c3_kernel!(
    gross_assimilation::AbstractVector{T},
    net_assimilation::AbstractVector{T},
    water_limited_assimilation::AbstractVector{T},
    leaf_respiration::AbstractVector{T},
    potential_vcmax::AbstractVector{T},
    vcmax::AbstractVector{T},
    nitrogen_limitation::AbstractVector{T},
    lambda::AbstractVector{T},
    temperature_stress::AbstractVector{T},
    apar::AbstractVector{T},
    daylength::AbstractVector{T},
    temperature::AbstractVector{T},
    co2::AbstractVector{T},
    fixed_parameters,
    comp_vcmax::Bool,
) where {T <: AbstractFloat}
    cell = @index(Global)
    parameters = kernel_value(fixed_parameters)
    CFT = parameters.CFT
    lpjmlparams = parameters.lpjmlparams
    photoparams = parameters.photoparams
    @unpack b = CFT
    @unpack ko25, kc25, alphac3, theta, LAMBDA_OPT = lpjmlparams
    @unpack q10ko, q10kc, po2, tau25, q10tau, cmass, cq, p, lambdamc3 = photoparams

    stress = temperature_stress[cell]
    inactive = stress < T(1e-2)
    temperature_cell = temperature[cell]
    co2_cell = co2[length(co2) == 1 ? 1 : cell]
    ko = T(ko25) * T(q10ko)^((temperature_cell - T(25)) * T(0.1))
    kc = T(kc25) * T(q10kc)^((temperature_cell - T(25)) * T(0.1))
    fac = kc * (one(T) + T(po2) / ko)
    tau = T(tau25) * T(q10tau)^((temperature_cell - T(25)) * T(0.1))
    gammastar = T(po2) / (T(2) * tau)

    if comp_vcmax
        internal_co2 = T(lambdamc3) * co2_cell
        c1 = stress * T(alphac3) *
            ((internal_co2 - gammastar) / (internal_co2 + T(2) * gammastar))
        c2 = (internal_co2 - gammastar) / (internal_co2 + fac)
        s = T(24) / daylength[cell] * T(b)
        sigma = one(T) - (c2 - s) / (c2 - T(theta) * s)
        sigma = sqrt(max(zero(T), sigma))
        lambda[cell] = T(LAMBDA_OPT)
        potential = (one(T) / T(b)) * (c1 / c2) *
            ((T(2) * T(theta) - one(T)) * s -
             (T(2) * T(theta) * s - c2) * sigma) *
            apar[cell] * T(cmass) * T(cq)
        vcmax[cell] = inactive ? zero(T) : max(zero(T), potential)
        potential_vcmax[cell] = vcmax[cell]
        nitrogen_limitation[cell] = vcmax[cell] > zero(T) ? one(T) : zero(T)
    end

    internal_co2 = lambda[cell] * co2_cell
    c1 = stress * T(alphac3) *
        ((internal_co2 - gammastar) / (internal_co2 + T(2) * gammastar))
    c2 = (internal_co2 - gammastar) / (internal_co2 + fac)
    je = c1 * apar[cell] * T(cmass) * T(cq) / (daylength[cell] + T(1e-5))
    jc = c2 * hour2day(vcmax[cell])
    agd = compute_co_limited_assimilation(je, jc, T(theta), daylength[cell])
    gross = inactive ? zero(T) : max(zero(T), agd)
    gross_assimilation[cell] = gross
    leaf = inactive ? zero(T) : T(b) * vcmax[cell]
    leaf_respiration[cell] = leaf
    net_assimilation[cell], adt = compute_net_assimilation(gross, leaf, daylength[cell])
    water_limited_assimilation[cell] = compute_water_limited_assimilation(
        adt, T(cmass), temperature_cell, T(p),
    )
end

"""Cell-local C4 photosynthesis kernel with no intermediate device arrays."""
function photosynthesis_C4!(CFT::CFTParameters,
                            crop,
                            apar::AbstractArray{T},
                            pet_daylength::AbstractArray{T},
                            temp::AbstractArray{T};
                            lpjmlparams::LPJmLParams = lpjmlparams,
                            photoparams::PhotoParams = photoparams,
                            comp_vcmax = false
) where {T <: AbstractFloat}
    launch_1D!(
        photosynthesis_c4_kernel!,
        crop_fluxes(crop).carbon.gross_assimilation,
        crop_fluxes(crop).carbon.net_assimilation,
        crop_fluxes(crop).carbon.water_limited_assimilation,
        crop_fluxes(crop).carbon.leaf_respiration,
        crop_photosynthesis_auxiliary(crop).potential_vcmax,
        crop_photosynthesis_auxiliary(crop).vcmax,
        crop_photosynthesis_auxiliary(crop).nitrogen_limitation,
        crop_photosynthesis_auxiliary(crop).lambda,
        crop_photosynthesis_auxiliary(crop).temperature_stress,
        apar,
        pet_daylength,
        temp,
        kernel_constant(
            crop_fluxes(crop).carbon.gross_assimilation,
            (; CFT, lpjmlparams, photoparams),
        ),
        comp_vcmax,
    )
    return nothing
end

"""Dispatch photosynthesis to the compile-time C3 or C4 pathway."""
photosynthesis!(::Val{:C3}, CFT, crop, apar, daylength, temperature, co2; kwargs...) =
    photosynthesis_C3!(CFT, crop, apar, daylength, temperature, co2; kwargs...)
photosynthesis!(::Val{:C4}, CFT, crop, apar, daylength, temperature, co2; kwargs...) =
    photosynthesis_C4!(CFT, crop, apar, daylength, temperature; kwargs...)

@kernel inbounds = true function photosynthesis_c4_kernel!(
    gross_assimilation::AbstractVector{T},
    net_assimilation::AbstractVector{T},
    water_limited_assimilation::AbstractVector{T},
    leaf_respiration::AbstractVector{T},
    potential_vcmax::AbstractVector{T},
    vcmax::AbstractVector{T},
    nitrogen_limitation::AbstractVector{T},
    lambda::AbstractVector{T},
    temperature_stress::AbstractVector{T},
    apar::AbstractVector{T},
    daylength::AbstractVector{T},
    temperature::AbstractVector{T},
    fixed_parameters,
    comp_vcmax::Bool,
) where {T <: AbstractFloat}
    cell = @index(Global)
    parameters = kernel_value(fixed_parameters)
    CFT = parameters.CFT
    lpjmlparams = parameters.lpjmlparams
    photoparams = parameters.photoparams
    @unpack b = CFT
    @unpack alphac4, theta, LAMBDA_OPT = lpjmlparams
    @unpack lambdamc4, cmass, cq, p = photoparams

    stress = temperature_stress[cell]
    inactive = stress < T(1e-2)
    if comp_vcmax
        c1 = stress * T(alphac4)
        s = T(24) / daylength[cell] * T(b)
        sigma = one(T) - (one(T) - s) / (one(T) - T(theta) * s)
        sigma = sqrt(max(zero(T), sigma))
        lambda[cell] = T(LAMBDA_OPT)
        potential = (one(T) / T(b)) * c1 *
            ((T(2) * T(theta) - one(T)) * s -
             (T(2) * T(theta) * s - one(T)) * sigma) *
            apar[cell] * T(cmass) * T(cq)
        vcmax[cell] = inactive ? zero(T) : max(zero(T), potential)
        potential_vcmax[cell] = vcmax[cell]
        nitrogen_limitation[cell] = vcmax[cell] > zero(T) ? one(T) : zero(T)
    end

    phipi = min(one(T), lambda[cell] / T(lambdamc4))
    c1 = stress * phipi * T(alphac4)
    je = c1 * apar[cell] * T(cmass) * T(cq) / (daylength[cell] + T(1e-5))
    jc = hour2day(vcmax[cell])
    agd = compute_co_limited_assimilation(je, jc, T(theta), daylength[cell])
    gross = inactive ? zero(T) : max(zero(T), agd)
    gross_assimilation[cell] = gross
    leaf = inactive ? zero(T) : T(b) * vcmax[cell]
    leaf_respiration[cell] = leaf
    net_assimilation[cell], adt = compute_net_assimilation(gross, leaf, daylength[cell])
    water_limited_assimilation[cell] = compute_water_limited_assimilation(
        adt, T(cmass), temperature[cell], T(p),
    )
end
