"""
petpar!(pet, day, lat, temp, lwnet, swdown; dayseconds=86400)

Compute daylength, PAR, and equilibrium evapotranspiration diagnostics.
"""

"""Allocation-free radiation and equilibrium-evaporation preprocessing."""
function petpar!(pet::PetPar,
                 day::Int64,
                 lat::AbstractArray{T},
                 temp::AbstractArray{T},
                 lwnet::AbstractArray{T},
                 swdown::AbstractArray{T};
                 dayseconds::Integer = 86400
) where {T <: AbstractFloat}
    launch_1D!(
        petpar_kernel!,
        pet.daylength,
        pet.par,
        pet.eeq,
        pet.albedo,
        lat,
        temp,
        lwnet,
        swdown,
        day,
        dayseconds,
    )
    return nothing
end

"""
    compute_daylength(latitude, solar_declination)

Return astronomical daylength in hours for one cell. Latitude deliberately
passes through Float64 trigonometry before the result narrows to model
precision, preserving the existing CPU/GPU-compatible numerical contract.
"""
@inline function compute_daylength(latitude::T,
                                   solar_declination::T) where {T <: AbstractFloat}
    latitude_radians = Float64(latitude) * π / 180.0
    u = T(sin(latitude_radians) * sin(solar_declination))
    v = T(cos(latitude_radians) * cos(solar_declination))
    return u >= v ? T(24) :
           u <= -v ? zero(T) :
           T(T(24) * acos(-u / v) * (1 / π))
end

"""
    compute_equilibrium_evaporation(temperature, shortwave, longwave, albedo,
                                    daylength, dayseconds)

Compute LPJmL's daily equilibrium evaporation diagnostic in mm day⁻¹.
LPJmL clips negative energy-balance values to zero, but does not impose an
arbitrary upper bound.
"""
@inline function compute_equilibrium_evaporation(temperature::T,
                                                  shortwave::T,
                                                  longwave::T,
                                                  albedo::T,
                                                  daylength::T,
                                                  dayseconds::T) where {T <: AbstractFloat}
    temperature_offset = T(237.3) + temperature
    slope = T(2.503e6) * exp(T(17.269) * temperature / temperature_offset) /
        (temperature_offset^2)
    psychrometric = T(65.05) + T(0.064) * temperature
    latent_heat = T(2.495e6) - T(2380) * temperature
    net_shortwave = (one(T) - albedo) * shortwave
    equilibrium = dayseconds * slope / (slope + psychrometric) / latent_heat *
        (net_shortwave + longwave * daylength / T(24))
    return max(equilibrium, zero(T))
end

"""
    compute_absorbed_fraction(actual_lai, light_extinction, maize, snow_free)

Compute green-canopy fPAR using the ordinary Beer--Lambert branch or LPJmL's
maize branch, then suppress it under snow cover.
"""
@inline function compute_absorbed_fraction(actual_lai::T,
                                           light_extinction::T,
                                           maize::Bool,
                                           snow_free::Bool) where {T <: AbstractFloat}
    fraction = maize ?
        min(one(T), max(zero(T), T(0.2558) * max(T(0.01), actual_lai) - T(0.0024))) :
        one(T) - exp(-light_extinction * actual_lai)
    return fraction * T(snow_free)
end

"""
    compute_absorbed_par(par, leaf_albedo, alpha_a, absorbed_fraction)

Convert absorbed radiation fraction to canopy APAR while keeping the
photosynthetic conversion factor explicit.
"""
@inline function compute_absorbed_par(par::T,
                                      leaf_albedo::T,
                                      alpha_a::T,
                                      absorbed_fraction::T) where {T <: AbstractFloat}
    return par * (one(T) - leaf_albedo) * alpha_a * absorbed_fraction
end

@kernel inbounds = true function petpar_kernel!(
    daylength::AbstractVector{T},
    par::AbstractVector{T},
    eeq::AbstractVector{T},
    albedo::AbstractVector{T},
    latitude::AbstractVector{T},
    temperature::AbstractVector{T},
    lwnet::AbstractVector{T},
    swdown::AbstractVector{T},
    day::Integer,
    dayseconds::Integer,
) where {T <: AbstractFloat}
    cell = @index(Global)
    delta = T(deg2rad(-23.4 * cos(2 * π * (day + 10) / 365)))
    seconds = T(dayseconds)
    daylight = compute_daylength(latitude[cell], delta)
    daylength[cell] = daylight
    par[cell] = seconds * swdown[cell] / T(2)
    eeq[cell] = compute_equilibrium_evaporation(
        temperature[cell], swdown[cell], lwnet[cell], albedo[cell], daylight, seconds,
    )
end

@kernel inbounds = true function daylength_kernel!(
                                   pet_daylength::AbstractArray{T},
                                   u::AbstractArray{T},
                                   v::AbstractArray{T}
) where {T <: AbstractFloat}

    cell = @index(Global)

    if u[cell] >= v[cell]
        pet_daylength[cell] = 24
    elseif u[cell] <= -v[cell]
        pet_daylength[cell] = 0
    else
        hh = acos(-u[cell] / v[cell])
        pet_daylength[cell] = 24 * hh * (1 / π)
    end
end

# for one cft
"""
apar_crop!(CFT, crop, pet)

Compute absorbed PAR and fPAR. Set `maize=true` for the maize-specific fPAR
parameterization.
"""

function apar_crop!(
    CFT::CFTParameters, crop, pet::PetPar, snow_height = nothing; maize::Bool = false,
)
    launch_1D!(
        apar_crop_kernel!,
        crop_canopy_auxiliary(crop).apar,
        crop_canopy_auxiliary(crop).fpar,
        crop_prognostic(crop).canopy.lai,
        crop_prognostic(crop).canopy.lai_npp_deficit,
        snow_height === nothing ? pet.eeq : snow_height,
        pet.par,
        kernel_constant(crop_canopy_auxiliary(crop).apar, CFT),
        maize,
        snow_height !== nothing,
    )
    return nothing
end
# Radiation and daylength preprocessing for canopy photosynthesis.

"""
apar_crop_maize!(CFT, crop, pet)

Compute absorbed PAR and maize-specific fPAR parameterization.
"""


apar_crop_maize!(CFT::CFTParameters, crop, pet::PetPar, snow_height = nothing) =
    apar_crop!(CFT, crop, pet, snow_height; maize = true)

@kernel inbounds = true function apar_crop_kernel!(
    apar::AbstractVector{T},
    fpar::AbstractVector{T},
    lai::AbstractVector{T},
    lai_npp_deficit::AbstractVector{T},
    snow_height::AbstractVector{T},
    par::AbstractVector{T},
    fixed_cft,
    maize::Bool,
    apply_snow_cover::Bool,
) where {T <: AbstractFloat}
    cell = @index(Global)
    CFT = kernel_value(fixed_cft)
    light_extinction = T(CFT.lightextcoeff)
    leaf_albedo = T(CFT.albedo_leaf)
    alpha_a = T(CFT.alphaa)
    actual_lai = max(zero(T), lai[cell] - lai_npp_deficit[cell])
    snow_free = !apply_snow_cover || snow_height[cell] <= zero(T)
    absorbed_fraction = compute_absorbed_fraction(
        actual_lai, light_extinction, maize, snow_free,
    )
    fpar[cell] = absorbed_fraction
    apar[cell] = compute_absorbed_par(
        par[cell], leaf_albedo, alpha_a, absorbed_fraction,
    )
end
