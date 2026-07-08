"""
carbon_allocation!(PFT, crop, photos)

Partition crop biomass among leaf/root/storage/pool carbon compartments.
"""
function carbon_allocation!(PFT::PftParameters,
                            crop::Crop,
                            photos::Photos
)
    # 1D cell-wise allocation; crop.stoc provides launch length and kernel arg #1.
    kernel_params = (FROOTMAX = 0.4f0, FROOTMIN = 0.3f0)

    launch_1D!(carbon_allocation_kernel!,
               crop.stoc,
               crop.isgrowing,
               crop.growingdays,
               crop.vscal_sum,
               crop.vscal,
               crop.ndf,
               crop.wdf,
               crop.fphu,
               crop.senescence,
               crop.biomass,
               crop.resp,
               photos.agd,
               photos.rd,
               crop.npp,
               crop.lai,
               crop.leafc,
               crop.rootc,
               crop.poolc,
               crop.lai_nppdeficit,
               PFT,
               kernel_params)

end

@kernel inbounds = true function carbon_allocation_kernel!(
                                           crop_stoc::AbstractArray{T},
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
                                           crop_leafc::AbstractArray{T},
                                           crop_rootc::AbstractArray{T},
                                           crop_poolc::AbstractArray{T},
                                           crop_lai_nppdeficit::AbstractArray{T},
                                           PFT::PftParameters,
                                           kernel_params
) where {T <: AbstractFloat, B <: Bool, S <: Integer}

    cell = @index(Global)

    @unpack sla, hiopt, himin = PFT
    @unpack FROOTMAX, FROOTMIN = kernel_params

    if crop_isgrowing[cell] == 1
        # Undo LAI deficit correction from previous step before current-day allocation.
        crop_lai[cell] = crop_lai[cell] - crop_lai_nppdeficit[cell]
        # NPP = gross daytime photosynthesis - dark respiration - growth respiration bookkeeping.
        crop_npp[cell] = (photos_agd[cell] - photos_rd[cell] - crop_resp[cell])
        if ((crop_biomass[cell] + crop_npp[cell]) <= T(0.0001)) || ((crop_lai[cell] <= zero(T)) && (!crop_senescence[cell]))
            crop_poolc[cell] += crop_npp[cell]
            crop_biomass[cell] += crop_npp[cell]
        else
            crop_biomass[cell] += crop_npp[cell]
            crop_vscal_sum[cell] += crop_vscal[cell]
            if crop_growingdays[cell] > 0
                crop_ndf[cell] = crop_vscal_sum[cell] / crop_growingdays[cell] * 100
            else
                crop_ndf[cell] = T(100)
            end

            # Root carbon follows SWAT-style stress-scaled partitioning.
            df = min(crop_wdf[cell], crop_ndf[cell])
            froot = FROOTMAX - (FROOTMIN * crop_fphu[cell]) * df / (df + exp(T(6.13) - T(0.0883) * df))
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
            fhiopt = 100 * crop_fphu[cell] / (100 * crop_fphu[cell] + exp(T(11.1) - T(10.0) * crop_fphu[cell]))
            hi = hiopt > 1.0 ? fhiopt * (hiopt - one(T)) + one(T) : fhiopt * hiopt
            himind = himin > 1.0 ? fhiopt * (himin - one(T)) + one(T) : fhiopt * himin

            if crop_wdf[cell] >= zero(T)
                hi = (hi - himind) * crop_wdf[cell] / (crop_wdf[cell] + exp(T(6.13) -T(0.0883) * crop_wdf[cell])) + himind
            end

            if (crop_leafc[cell] + crop_rootc[cell]) < crop_biomass[cell]
                if hiopt > 1.0
                    crop_stoc[cell] = (one(T) - one(T) / hi) * (one(T) - froot) * crop_biomass[cell]
                else
                    crop_stoc[cell] = hi * (one(T) - froot) * crop_biomass[cell]
                end
                if (crop_leafc[cell] + crop_rootc[cell] + crop_stoc[cell]) > crop_biomass[cell]
                    crop_stoc[cell] = crop_biomass[cell] - crop_leafc[cell] - crop_rootc[cell]
                end
            else
                crop_stoc[cell] = zero(T)
            end

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
    end

end
