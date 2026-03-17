function transpiration!(photos_adtmm::AbstractArray{T},
                        param::LPJmLParam,
                        PFT::PftParameters,
                        crop::Crop,
                        pet::PetPar,
                        soil::Soil,
                        co2::AbstractArray{T}
) where {T <: AbstractFloat}

    @unpack LAMBDA_OPT = param

    crop.gp = (1.6f0 * photos_adtmm ./ (ppm2bar(co2) * (1 - LAMBDA_OPT) .* hour2sec(pet.daylength))) .+ crop.fpar # potential canopy conductance

    wr = sum(soil.w .* crop.rootdist, dims = 1)
    # supply = emax * wr .* (1 .- exp.(-0.04f0 * crop.rootc))
    # demand = ifelse.(crop.gp .> 0, (1 .- crop.canopy_wet) .* pet.eeq * ALPHAM ./ (1 .+ (GM * ALPHAM) ./ crop.gp), zero(T))
    # transp = ifelse.(wr .> 0, min.(supply, demand) ./ wr .* fpc, zero(T)) # here the crop.fpc = 1, so we just omit it in the kernel fucntion

    backend = get_backend(crop.gp)

    kernel = water_demand_supply_kernel!(backend)
    
    kernel(param, 
           PFT, 
           crop.trans_layer, 
           crop.w_demandsum,
           crop.w_supplysum,
           crop.wdf,
           crop.wscal,
           crop.gp, 
           crop.rootc, 
           crop.canopy_wet,
           crop.isgrowing,
           pet.eeq, 
           crop.rootdist, 
           soil.w, 
           soil.whcs, 
           wr, 
           ndrange=length(crop.gp))

    synchronize(backend)

end

@kernel function water_demand_supply_kernel!(param::LPJmLParam,
                                             PFT::PftParameters,
                                             crop_trans_layer::AbstractArray{T},
                                             crop_w_demandsum::AbstractArray{T},
                                             crop_w_supplysum::AbstractArray{T},
                                             crop_wdf::AbstractArray{T},
                                             crop_wscal::AbstractArray{T},
                                             crop_gp::AbstractArray{T},
                                             crop_rootc::AbstractArray{T},
                                             crop_canopy_wet::AbstractArray{T},
                                             crop_isgrowing::AbstractArray{S},
                                             pet_eeq::AbstractArray{T},
                                             crop_rootdist::AbstractArray{T},
                                             soil_w::AbstractArray{M},
                                             soil_whcs::AbstractArray{M},
                                             wr::AbstractArray{T};
                                             soil_layers = 5
) where {T <: AbstractFloat, M <: AbstractFloat, S <: Integer}
    
    cell = @index(Global)

    @unpack ALPHAM, GM = param
    @unpack fpc, emax = PFT

    if crop_isgrowing[cell] == 1
        supply = emax * wr[cell] * (1 - exp(-0.04f0 * crop_rootc[cell]))
        if crop_gp[cell] > 0
            demand = (1 - crop_canopy_wet[cell]) * pet_eeq[cell] * ALPHAM / (1 + (GM * ALPHAM) / crop_gp[cell])
        else
            demand = zero(T)
        end
        
        crop_w_demandsum[cell] += demand
        if supply > demand
            crop_w_supplysum[cell] += demand
        else
            crop_w_supplysum[cell] += supply
        end

        if crop_w_demandsum[cell] > 0.0
            crop_wdf[cell] = T(100.0) * crop_w_supplysum[cell] / crop_w_demandsum[cell]
        else
            crop_wdf[cell] = T(100.0)
        end

        if pet_eeq[cell] > 0.0 && crop_gp[cell] > 0.0
            crop_wscal[cell] = (emax * wr[cell]) / (pet_eeq[cell] * ALPHAM / (one(T) + (GM * ALPHAM) / crop_gp[cell]))
            if crop_wscal[cell] > 1.0
                crop_wscal[cell] = one(T)
            end
        else
            crop_wscal[cell] = one(T)
        end

        if wr[cell] > 0
            transp = min(supply, demand) / wr[cell] * fpc
        else
            transp = zero(T)
        end

        transp_cor = zero(T)

        if transp > 0
            for l in 1:soil_layers
                transp_frac = 1
                if transp * crop_rootdist[l] * soil_w[l, cell] > soil_w[l, cell] * soil_whcs[l, cell]
                    transp_frac = soil_whcs[l, cell] / (transp * crop_rootdist[l])
                end
                transp_tmp = transp * crop_rootdist[l] * soil_w[l, cell] * transp_frac
                if transp_tmp > soil_w[l, cell] * soil_whcs[l, cell]
                    transp_cor += soil_w[l, cell] * soil_whcs[l, cell]
                    if transp_cor < 1.0f-5
                        transp_cor = zero(T)
                    end
                else
                    transp_cor += transp_tmp
                end
            end
        else
            transp_cor = zero(T)
        end

        if wr[cell] > 0
            transp = transp_cor / wr[cell]
        else
            transp = zero(T)
        end

        for l in 1:soil_layers
            crop_trans_layer[l, cell] = transp * crop_rootdist[l] * soil_w[l, cell]
            if crop_trans_layer[l, cell] > soil_w[l, cell] * soil_whcs[l, cell]
                crop_trans_layer[l, cell] = soil_w[l, cell] * soil_whcs[l, cell]
            end
        end
    else
        for l in 1:soil_layers
            crop_trans_layer[l, cell] = zero(T)
        end
        crop_w_demandsum[cell] = zero(T)
        crop_w_supplysum[cell] = zero(T)
        crop_wdf[cell] = zero(T)
        crop_wscal[cell] = zero(T)
    end
end