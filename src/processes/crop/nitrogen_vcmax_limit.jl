"""
    limit_vcmax_by_nitrogen!(crop, CFT, temperature)

Apply LPJmL's crop leaf-nitrogen constraint to the potential Rubisco capacity.
`crop_stress_auxiliary(crop).nitrogen_demand_leaf` is the leaf-N stock remaining after uptake logic;
structural leaf N is protected and only excess N supports Rubisco activity.
The result is bounded to `[eps(T), potential_vcmax]` for an active crop and the
dimensionless retained fraction is stored in
`crop_photosynthesis_auxiliary(crop).nitrogen_limitation`.
"""
function limit_vcmax_by_nitrogen!(crop,
                                 CFT::CFTParameters,
                                 temperature::AbstractArray{T};
                                 lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat}
    launch_1D!(
        nitrogen_vcmax_limit_kernel!,
        crop_photosynthesis_auxiliary(crop).vcmax,
        crop_photosynthesis_auxiliary(crop).potential_vcmax,
        crop_photosynthesis_auxiliary(crop).nitrogen_limitation,
        crop_stress_auxiliary(crop).nitrogen_demand_leaf,
        crop_prognostic(crop).carbon.leaf,
        crop_prognostic(crop).phenology.is_growing,
        temperature,
        CFT,
        lpjmlparams,
    )
    return nothing
end

@kernel inbounds = true function nitrogen_vcmax_limit_kernel!(
    vcmax::AbstractArray{T},
    potential_vcmax::AbstractArray{T},
    nitrogen_limitation::AbstractArray{T},
    available_leaf_nitrogen::AbstractArray{T},
    leaf_carbon::AbstractArray{T},
    is_growing::AbstractArray{S},
    temperature::AbstractArray{T},
    CFT::CFTParameters,
    lpjmlparams::LPJmLParams,
) where {T <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    potential = max(zero(T), potential_vcmax[cell])

    if is_growing[cell] == one(S) && potential > zero(T)
        @unpack ncleaf = CFT
        @unpack p, k_temp = lpjmlparams
        rubisco_nitrogen = max(
            zero(T),
            available_leaf_nitrogen[cell] - T(ncleaf.low) * leaf_carbon[cell],
        )
        nitrogen_capacity = rubisco_nitrogen /
            exp(-T(k_temp) * (temperature[cell] - T(25))) /
            (T(p) * T(1e-3)) * (T(86400) * T(12) * T(1e-6))
        limited = min(potential, max(eps(T), nitrogen_capacity))
        vcmax[cell] = limited
        nitrogen_limitation[cell] = clamp(limited / potential, zero(T), one(T))
    else
        vcmax[cell] = potential
        nitrogen_limitation[cell] = potential > zero(T) ? one(T) : zero(T)
    end
end
