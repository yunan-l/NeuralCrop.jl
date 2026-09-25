"""
carbon_allocation!(CFT, crop, photos)

Partition crop biomass among leaf/root/storage/pool carbon compartments.
"""
function carbon_allocation!(CFT::CFTParameters,
                            crop
)
    # 1D cell-wise allocation; crop_prognostic(crop).carbon.storage provides launch length and kernel arg #1.
    T = eltype(crop_prognostic(crop).carbon.storage)
    kernel_params = (FROOTMAX = T(0.4), FROOTMIN = T(0.3))

    launch_1D!(carbon_allocation_kernel!,
               crop_prognostic(crop).carbon.storage,
               crop_events(crop).harvest,
               crop_prognostic(crop).phenology.is_growing,
               crop_prognostic(crop).phenology.growing_days,
               crop_prognostic(crop).nitrogen.stress_sum,
               crop_prognostic(crop).nitrogen.sufficiency,
               crop_stress_auxiliary(crop).nitrogen_deficit,
               crop_stress_auxiliary(crop).water_deficit,
               crop_phenology_auxiliary(crop).fphu,
               crop_prognostic(crop).phenology.senescence,
               crop_prognostic(crop).carbon.biomass,
               crop_fluxes(crop).carbon.respiration,
               crop_fluxes(crop).carbon.gross_assimilation,
               crop_fluxes(crop).carbon.leaf_respiration,
               crop_fluxes(crop).carbon.npp,
               crop_prognostic(crop).canopy.lai,
               crop_canopy_auxiliary(crop).actual_lai,
               crop_prognostic(crop).carbon.leaf,
               crop_prognostic(crop).carbon.root,
               crop_prognostic(crop).carbon.pool,
               crop_prognostic(crop).canopy.lai_npp_deficit,
               kernel_constant(
                   crop_prognostic(crop).carbon.storage,
                   (; CFT, kernel_params),
               ))

end

"""Compute daily NPP after leaf, maintenance, and growth respiration."""
@inline compute_crop_npp(
    gross_assimilation::T, leaf_respiration::T, crop_respiration::T,
) where {T <: AbstractFloat} =
    gross_assimilation - leaf_respiration - crop_respiration

"""Compute the seasonal nitrogen sufficiency percentage used by root allocation."""
@inline function compute_seasonal_nitrogen_sufficiency(
    accumulated_sufficiency::T, growing_days::S,
) where {T <: AbstractFloat, S <: Integer}
    return growing_days > zero(S) ? accumulated_sufficiency / T(growing_days) * T(100) : T(100)
end

"""Compute SWAT-style root-carbon fraction from water/N stress and phenology."""
@inline function compute_root_carbon_fraction(
    water_sufficiency::T,
    nitrogen_sufficiency::T,
    phenology_fraction::T,
    root_maximum::T,
    root_minimum::T,
) where {T <: AbstractFloat}
    stress = min(water_sufficiency, nitrogen_sufficiency)
    return root_maximum - (root_minimum * phenology_fraction) * stress /
           (stress + exp(T(6.13) - T(0.0883) * stress))
end

"""Compute water-limited LPJmL harvest index for one crop stand."""
@inline function compute_harvest_index(
    phenology_fraction::T,
    optimal_index::T,
    minimum_index::T,
    water_sufficiency::T,
) where {T <: AbstractFloat}
    potential = T(100) * phenology_fraction /
                (T(100) * phenology_fraction +
                 exp(T(11.1) - T(10) * phenology_fraction))
    optimal = optimal_index > one(T) ? potential * (optimal_index - one(T)) + one(T) :
              potential * optimal_index
    minimum = minimum_index > one(T) ? potential * (minimum_index - one(T)) + one(T) :
              potential * minimum_index
    water_sufficiency >= zero(T) || return optimal
    return (optimal - minimum) * water_sufficiency /
           (water_sufficiency + exp(T(6.13) - T(0.0883) * water_sufficiency)) + minimum
end

"""Compute and mass-cap storage carbon after leaf/root allocation."""
@inline function compute_storage_carbon(
    biomass::T,
    leaf_carbon::T,
    root_carbon::T,
    root_fraction::T,
    harvest_index::T,
    optimal_index::T,
) where {T <: AbstractFloat}
    leaf_carbon + root_carbon < biomass || return zero(T)
    candidate = optimal_index > one(T) ?
                (one(T) - one(T) / harvest_index) * (one(T) - root_fraction) * biomass :
                harvest_index * (one(T) - root_fraction) * biomass
    return min(candidate, biomass - leaf_carbon - root_carbon)
end

@kernel inbounds = true function carbon_allocation_kernel!(
                                           crop_stoc::AbstractArray{T},
                                           crop_harvest::AbstractArray{S},
                                           crop_isgrowing::AbstractArray{S},
                                           crop_growingdays::AbstractArray{S},
                                           crop_vscal_sum::AbstractArray{T},
                                           crop_vscal::AbstractArray{T},
                                           crop_ndf::AbstractArray{T},
                                           crop_wdf::AbstractArray{T},
                                           crop_fphu::AbstractArray{T},
                                           crop_senescence::AbstractArray{B},
                                           crop_biomass::AbstractArray{T},
                                           crop_resp::AbstractArray{T},
                                           photos_agd::AbstractArray{T},
                                           photos_rd::AbstractArray{T},
                                           crop_npp::AbstractArray{T},
                                           crop_lai::AbstractArray{T},
                                           crop_actual_lai::AbstractArray{T},
                                           crop_leafc::AbstractArray{T},
                                           crop_rootc::AbstractArray{T},
                                           crop_poolc::AbstractArray{T},
                                           crop_lai_nppdeficit::AbstractArray{T},
                                           fixed_parameters,
) where {T <: AbstractFloat, B <: Bool, S <: Integer}

    cell = @index(Global)
    parameters = kernel_value(fixed_parameters)
    CFT = parameters.CFT
    kernel_params = parameters.kernel_params

    @unpack sla, hiopt, himin = CFT
    @unpack FROOTMAX, FROOTMIN = kernel_params

    if crop_isgrowing[cell] == 1
        # LPJmL preserves the potential phenological LAI and applies the NPP
        # deficit only when actual LAI is consumed or reported.
        actual_lai = max(zero(T), crop_lai[cell] - crop_lai_nppdeficit[cell])
        # Complete crop carbon cost: leaf respiration plus maintenance/growth
        # respiration, including root respiration.
        crop_npp[cell] = compute_crop_npp(photos_agd[cell], photos_rd[cell], crop_resp[cell])
        if ((crop_biomass[cell] + crop_npp[cell]) <= T(0.0001)) || ((actual_lai <= zero(T)) && (!crop_senescence[cell]))
            # LPJmL reports `negbm` here. The daily driver then harvests the
            # remaining pools and removes the failed crop stand.
            crop_poolc[cell] += crop_npp[cell]
            crop_biomass[cell] += crop_npp[cell]
            crop_harvest[cell] = one(S)
        else
            crop_biomass[cell] += crop_npp[cell]
            crop_vscal_sum[cell] += crop_vscal[cell]
            crop_ndf[cell] = compute_seasonal_nitrogen_sufficiency(
                crop_vscal_sum[cell], crop_growingdays[cell],
            )

            # Root carbon follows SWAT-style stress-scaled partitioning.
            froot = compute_root_carbon_fraction(
                crop_wdf[cell], crop_ndf[cell], crop_fphu[cell], FROOTMAX, FROOTMIN,
            )
            crop_rootc[cell] = froot * crop_biomass[cell]

            # Leaf carbon is constrained by LAI and SLA; in senescence it is mass-balanced.
            if !crop_senescence[cell]
                if (crop_biomass[cell] - crop_rootc[cell]) >= (crop_lai[cell] / sla)
                    crop_leafc[cell] = crop_lai[cell] / sla
                    crop_lai_nppdeficit[cell] = zero(T)
                else
                    crop_leafc[cell] = crop_biomass[cell] - crop_rootc[cell]
                    crop_lai_nppdeficit[cell] = crop_lai[cell] - crop_leafc[cell] * sla
                end
            else
                if (crop_leafc[cell] + crop_rootc[cell] + crop_stoc[cell]) > crop_biomass[cell]
                    crop_leafc[cell] = crop_biomass[cell] - crop_rootc[cell] - crop_stoc[cell]
                end
                if crop_leafc[cell] < zero(T)
                    crop_leafc[cell] = zero(T)
                end
            end

            # Storage carbon (harvest index branch) is computed after leaf/root partitioning.
            hi = compute_harvest_index(crop_fphu[cell], T(hiopt), T(himin), crop_wdf[cell])
            crop_stoc[cell] = compute_storage_carbon(
                crop_biomass[cell], crop_leafc[cell], crop_rootc[cell], froot, hi, T(hiopt),
            )

            # Pool carbon closes biomass balance and is clipped during senescence if negative.
            crop_poolc[cell] = crop_biomass[cell] - crop_leafc[cell] - crop_rootc[cell] - crop_stoc[cell]
            # pool can become negative during senescence
            if crop_senescence[cell] && crop_poolc[cell] < zero(T)
                if (crop_stoc[cell] + crop_poolc[cell]) < zero(T)
                    crop_poolc[cell] += crop_stoc[cell]
                    crop_stoc[cell] = zero(T)
                    if (crop_rootc[cell] + crop_poolc[cell]) < zero(T)
                        crop_poolc[cell] += crop_rootc[cell]
                        crop_rootc[cell] = zero(T) # remainder negative pool must be compensated by leaves,
                        crop_leafc[cell] += crop_poolc[cell]
                        crop_poolc[cell] = zero(T)
                    else
                        crop_rootc[cell] += crop_poolc[cell]
                        crop_poolc[cell] = zero(T)
                    end
                else
                    crop_stoc[cell] += crop_poolc[cell]
                    crop_poolc[cell] = zero(T)
                end
            end
        end

    else
        crop_leafc[cell] = zero(T)
        crop_rootc[cell] = zero(T)
        crop_stoc[cell] = zero(T)
        crop_poolc[cell] = zero(T)
        crop_npp[cell] = zero(T)
        crop_biomass[cell] = zero(T)
        crop_vscal_sum[cell] = zero(T)
        crop_ndf[cell] = zero(T)
        crop_lai_nppdeficit[cell] = zero(T)
    end

    crop_actual_lai[cell] = max(zero(T), crop_lai[cell] - crop_lai_nppdeficit[cell])

end
