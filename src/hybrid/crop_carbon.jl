### crop carbon allocation with Neural ODE
function crop_carbon_node!(nn_model,
                           ps,
                           st,
                           photos::Photos,
                           crop::Crop,
                           PFT::PftParameters,
                           temp::AbstractArray{T},
                           temp_n::AbstractArray{T},
                           soil_swc::AbstractArray{M}
) where {T <: AbstractFloat, M <: AbstractFloat}

    # compute crop respiration
    Zygote.ignore() do
        respiration!(crop, PFT, temp, photos.agd - photos.rd)
    end

    # compute crop carbon allocation
    crop.npp = (photos.agd - photos.rd - crop.resp)
    crop.biomass = crop.biomass .+ crop.npp
    crop.biomass = crop.biomass .* crop.isgrowing
    
    # input = vcat(reshape(crop.npp/20, (1, :)), reshape(crop.lai/5, (1, :)), reshape(crop.leafn/5, (1, :)), mean(soil_swc[1:3, :], dims = 1), reshape(temp_n, (1, :))) .* reshape(crop.isgrowing, (1, :))
    input = vcat(reshape(crop.npp/20, (1, :)), reshape(crop.lai/5, (1, :)), reshape(crop.fphu, (1, :)), mean(soil_swc[1:3, :], dims = 1)) .* reshape(crop.isgrowing, (1, :))
    crop.vegc = neural_allocation(nn_model, crop.vegc, ps, st, input)
    
    crop.rootc = crop.vegc[1, :] .* crop.isgrowing
    crop.leafc = crop.vegc[2, :] .* crop.isgrowing
    crop.stoc = crop.vegc[3, :] .* crop.isgrowing
    crop.poolc = crop.vegc[4, :] .* crop.isgrowing
end


### crop carbon allocation with hybrid Neural ODE
function crop_carbon_hybrid!(nn_model,
                             ps,
                             st,
                             photos::Photos,
                             crop::Crop,
                             PFT::PftParameters,
                             temp::AbstractArray{T},
                             temp_n::AbstractArray{T},
                             soil_swc::AbstractArray{M};
                             lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat, M <: AbstractFloat}

    # compute crop respiration
    Zygote.ignore() do
        respiration!(crop, PFT, temp, photos.agd - photos.rd)
    end

    # compute crop root and leaf carbon allocation
    Zygote.ignore() do
        carbon_allocation_root_leaf!(PFT, crop, photos)
    end

    input = vcat(reshape(crop.npp/20, (1, :)), reshape(crop.fphu, (1, :)), reshape(temp_n, (1, :)), mean(soil_swc[1:3, :], dims = 1)) .* reshape(crop.isgrowing, (1, :))
    crop.stoc = neural_stoc(nn_model, reshape(crop.stoc/20, (1, :)), ps, st, input)

    # compute crop rest carbon allocation
    Zygote.ignore() do
        carbon_allocation_pool!(crop)
    end

    crop.vegc = vcat(reshape(crop.rootc, (1, :)), reshape(crop.leafc, (1, :)), reshape(crop.stoc, (1, :)), reshape(crop.poolc, (1, :)))

end

function carbon_allocation_root_leaf!(PFT::PftParameters,
                                      crop::Crop,
                                      photos::Photos
)

    backend = KernelAbstractions.get_backend(crop.stoc)

    kernel = carbon_allocation_leaf_root_kernel!(backend)
    
    kernel(PFT,
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
           crop.stoc,
           crop.poolc,
           crop.lai_nppdeficit;
           ndrange=length(crop.stoc))
    
    KernelAbstractions.synchronize(backend)

end

@kernel function carbon_allocation_leaf_root_kernel!(PFT::PftParameters,
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
                                                     crop_stoc::AbstractArray{T},
                                                     crop_poolc::AbstractArray{T},
                                                     crop_lai_nppdeficit::AbstractArray{T};
                                                     FROOTMAX = 0.4f0,
                                                     FROOTMIN = 0.3f0
) where {T <: AbstractFloat, B <: Bool, S <: Integer}

    cell = @index(Global)

    @unpack sla = PFT
    # @unpack sla, hiopt, himin = PFT

    if crop_isgrowing[cell] == 1
        crop_lai[cell] = crop_lai[cell] - crop_lai_nppdeficit[cell]
        crop_npp[cell] = (photos_agd[cell] - photos_rd[cell] - crop_resp[cell])
        if ((crop_biomass[cell] + crop_npp[cell]) <= 0.0001) || ((crop_lai[cell] <= 0.0) && (!crop_senescence[cell]))
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

            # root carbon
            df = min(crop_wdf[cell], crop_ndf[cell])
            froot = FROOTMAX - (FROOTMIN * crop_fphu[cell]) * df / (df + exp(T(6.13) - T(0.0883) * df))
            crop_rootc[cell] = froot * crop_biomass[cell]

            # leaf carbon
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
                if crop_leafc[cell] < 0
                    crop_leafc[cell] = zero(T)
                end
            end

            # # storage carbon
            # fhiopt = 100 * crop_fphu[cell] / (100 * crop_fphu[cell] + exp(T(11.1) - T(10.0) * crop_fphu[cell]))
            # hi = hiopt > 1.0 ? fhiopt * (hiopt - one(T)) + one(T) : fhiopt * hiopt
            # himind = himin > 1.0 ? fhiopt * (himin - one(T)) + one(T) : fhiopt * himin

            # if crop_wdf[cell] >= 0
            #     hi = (hi - himind) * crop_wdf[cell] / (crop_wdf[cell] + exp(T(6.13) -T(0.0883) * crop_wdf[cell])) + himind
            # end

            # if (crop_leafc[cell] + crop_rootc[cell]) < crop_biomass[cell]
            #     if hiopt > 1.0
            #         crop_stoc[cell] = (one(T) - one(T) / hi) * (one(T) - froot) * crop_biomass[cell]
            #     else
            #         crop_stoc[cell] = hi * (one(T) - froot) * crop_biomass[cell]
            #     end
            #     if (crop_leafc[cell] + crop_rootc[cell] + crop_stoc[cell]) > crop_biomass[cell]
            #         crop_stoc[cell] = crop_biomass[cell] - crop_leafc[cell] - crop_rootc[cell]
            #     end
            # else
            #     crop_stoc[cell] = zero(T)
            # end

            # # pool carbon
            # crop_poolc[cell] = crop_biomass[cell] - crop_leafc[cell] - crop_rootc[cell] - crop_stoc[cell]
            # # pool can become negative during senescence
            # if crop_senescence[cell] && crop_poolc[cell] < 0.0
            #     if (crop_stoc[cell] + crop_poolc[cell]) < 0.0
            #         crop_poolc[cell] += crop_stoc[cell]
            #         crop_stoc[cell] = zero(T)
            #         if (crop_rootc[cell] + crop_poolc[cell]) < 0.0
            #             crop_poolc[cell] += crop_rootc[cell]
            #             crop_rootc[cell] = zero(T) # remainder negative pool must be compensated by leaves,
            #             crop_leafc[cell] += crop_poolc[cell]
            #             crop_poolc[cell] = zero(T)
            #         else
            #             crop_rootc[cell] += crop_poolc[cell]
            #             crop_poolc[cell] = zero(T)
            #         end
            #     else
            #         crop_stoc[cell] += crop_poolc[cell]
            #         crop_poolc[cell] = zero(T)
            #     end
            # end
        end

    else
        crop_leafc[cell] = zero(T)
        crop_rootc[cell] = zero(T)
        # crop_poolc[cell] = zero(T)
        # crop_stoc[cell] = zero(T)
        crop_vscal_sum[cell] = zero(T)
        crop_ndf[cell] = zero(T)
    end
end


function carbon_allocation_pool!(crop::Crop)

    backend = KernelAbstractions.get_backend(crop.stoc)

    kernel = carbon_allocation_pool_kernel!(backend)
    
    kernel(crop.isgrowing,
           crop.senescence,
           crop.biomass,
           crop.leafc,
           crop.rootc,
           crop.stoc,
           crop.poolc;
           ndrange=length(crop.stoc))
    
    KernelAbstractions.synchronize(backend)
end

@kernel function carbon_allocation_pool_kernel!(crop_isgrowing::AbstractArray{S},
                                                crop_senescence::AbstractArray{B},
                                                crop_biomass::AbstractArray{T},
                                                crop_leafc::AbstractArray{T},
                                                crop_rootc::AbstractArray{T},
                                                crop_stoc::AbstractArray{T},
                                                crop_poolc::AbstractArray{T}
) where {T <: AbstractFloat, B <: Bool, S <: Integer}

    cell = @index(Global)

    if crop_isgrowing[cell] == 1
            # pool carbon
            crop_poolc[cell] = crop_biomass[cell] - crop_leafc[cell] - crop_rootc[cell] - crop_stoc[cell]
            # pool can become negative during senescence
            if crop_senescence[cell] && crop_poolc[cell] < 0.0
                if (crop_rootc[cell] + crop_poolc[cell]) < 0.0
                    crop_poolc[cell] += crop_rootc[cell]
                    crop_rootc[cell] = zero(T) # remainder negative pool must be compensated by leaves,
                    crop_leafc[cell] += crop_poolc[cell]
                    crop_poolc[cell] = zero(T)
                else
                    crop_rootc[cell] += crop_poolc[cell]
                    crop_poolc[cell] = zero(T)
                end
            end
    else
        crop_poolc[cell] = zero(T)
        crop_biomass[cell] = zero(T)
    end
end