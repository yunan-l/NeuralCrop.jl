"""
ndemand_crop!(crop, CFT, photos_vcmax, temp)

Compute crop nitrogen demand from photosynthetic potential and organ stoichiometry.
"""
function ndemand_crop!(crop,
                       CFT::CFTParameters,
                       photos_vcmax::AbstractArray{T},
                       temp::AbstractArray{T};
                       lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat}

    kernel_params = (lpjmlparams = lpjmlparams,)

    launch_1D!(
        ndemand_crop_kernel!,
        crop_stress_auxiliary(crop).nitrogen_demand_total,
        crop_prognostic(crop).carbon.leaf,
        crop_prognostic(crop).carbon.root,
        crop_prognostic(crop).carbon.pool,
        crop_prognostic(crop).carbon.storage,
        crop_stress_auxiliary(crop).nitrogen_demand_leaf,
        crop_prognostic(crop).phenology.is_growing,
        photos_vcmax,
        temp,
        kernel_constant(
            crop_stress_auxiliary(crop).nitrogen_demand_total,
            (; CFT, kernel_params),
        ),
    )

end

@kernel inbounds = true function ndemand_crop_kernel!(
                                      crop_ndemand_tot::AbstractArray{T},
                                      crop_leafc::AbstractArray{T},
                                      crop_rootc::AbstractArray{T},
                                      crop_poolc::AbstractArray{T},
                                      crop_stoc::AbstractArray{T},
                                      crop_ndemand_leaf::AbstractArray{T},
                                      crop_isgrowing::AbstractArray{S},
                                      photos_vcmax::AbstractArray{T},
                                      temp::AbstractArray{T},
                                      fixed_parameters,
) where {T <: AbstractFloat, S <: Integer}

    cell = @index(Global)
    parameters = kernel_value(fixed_parameters)
    CFT = parameters.CFT
    kernel_params = parameters.kernel_params

    @unpack lpjmlparams = kernel_params

    @unpack p, k_temp = lpjmlparams
    @unpack ratio, ncleaf = CFT

    if crop_isgrowing[cell] == 1
        # LPJmL ndemand_crop: Rubisco requirement plus structural minimum leaf N.
        rubisco_demand = T(p) * T(1e-3) * photos_vcmax[cell] /
                         (T(86400) * T(12) * T(1e-6)) *
                         exp(-T(k_temp) * (temp[cell] - T(25)))
        crop_ndemand_leaf[cell] = rubisco_demand + T(ncleaf.low) * crop_leafc[cell]

        nc_ratio = zero(T)
        if crop_leafc[cell] > zero(T)
            nc_ratio = crop_ndemand_leaf[cell] / crop_leafc[cell]
        end

        if nc_ratio > ncleaf.high
            nc_ratio = ncleaf.high
        elseif nc_ratio < ncleaf.low
            nc_ratio = ncleaf.low
        end

        crop_ndemand_tot[cell] = crop_ndemand_leaf[cell] + nc_ratio * (crop_rootc[cell] / ratio.root + crop_poolc[cell] / ratio.pool + crop_stoc[cell] / ratio.sto)
    else
        crop_ndemand_tot[cell] = zero(T)
        crop_ndemand_leaf[cell] = zero(T)
    end

end
