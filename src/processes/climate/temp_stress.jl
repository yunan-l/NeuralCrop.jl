"""
temp_stress(CFT, pet, photos, temp)

Compute temperature stress scalar used by photosynthesis routines.
"""
function temp_stress(CFT::CFTParameters,
                     pet::PetPar,
                     crop,
                     temp::AbstractArray{T};
                     photoparams::PhotoParams = photoparams
) where {T <: AbstractFloat}

    launch_1D!(
        temp_stress_kernel!,
        crop_photosynthesis_auxiliary(crop).temperature_stress,
        pet.daylength,
        temp,
        kernel_constant(
            crop_photosynthesis_auxiliary(crop).temperature_stress,
            (; CFT, photoparams),
        ),
    )

end

"""
    compute_photosynthesis_temperature_stress(daylength, temperature, path, ...)

LPJmL's smooth lower and upper temperature response for C3/C4 assimilation.
The function is scalar and allocation-free so it is safe to call from a CPU or
GPU kernel. `path` uses the existing `1 = C3`, `2 = C4` convention.
"""
@inline function compute_photosynthesis_temperature_stress(
    daylength::T,
    temperature::T,
    path,
    temperature_co2,
    temperature_photosynthesis,
    c3_maximum::T,
    c4_maximum::T,
) where {T <: AbstractFloat}
    # No light, or pathway-specific hard thermal limit: no canopy assimilation.
    (daylength < 0.01 ||
     (path == 1 && temperature > c3_maximum) ||
     (path == 2 && temperature > c4_maximum) ||
     temperature >= temperature_co2.high) && return zero(T)

    lower_slope = T(2 * log(1 / 0.99 - 1)) /
        (temperature_co2.low - temperature_photosynthesis.low)
    # LPJmL fscancftpar.c: midpoint of lower CO2 and photosynthesis limits.
    lower_midpoint = (T(temperature_co2.low) + T(temperature_photosynthesis.low)) * T(0.5)
    upper_slope = T(log(0.99 / 0.01)) /
        (temperature_co2.high - temperature_photosynthesis.high)
    lower_response = 1 / (1 + exp(lower_slope * (lower_midpoint - temperature)))
    upper_response = 1 - 0.01 * exp(
        upper_slope * (temperature - temperature_photosynthesis.high),
    )
    return T(lower_response * upper_response)
end


@kernel inbounds = true function temp_stress_kernel!(
                                     photos_tstress::AbstractArray{T},
                                     pet_daylength::AbstractArray{T},
                                     temp::AbstractArray{T},
                                     fixed_parameters,
) where {T <: AbstractFloat}

    cell = @index(Global)
    parameters = kernel_value(fixed_parameters)
    CFT = parameters.CFT
    photoparams = parameters.photoparams

    @unpack path, temp_co2, temp_photos = CFT
    @unpack tmc3, tmc4 = photoparams

    photos_tstress[cell] = compute_photosynthesis_temperature_stress(
        pet_daylength[cell], temp[cell], path, temp_co2, temp_photos,
        T(tmc3), T(tmc4),
    )
end
