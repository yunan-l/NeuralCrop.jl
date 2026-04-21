function nuptake_crop!(crop::Crop,
                       PFT::PftParameters,
                       soil::Soil;
                       lpjmlparams::LPJmLParams = lpjmlparams
)

    backend = KernelAbstractions.get_backend(crop.nitrogen)

    kernel = nuptake_crop_kernel!(backend)
    
    kernel(PFT,
           lpjmlparams,
           crop.leafn,
           crop.leafc,
           crop.rootn,
           crop.rootc,
           crop.ndemand_leaf,
           crop.ndemand_tot,
           crop.nitrogen,
           crop.vscal,
           crop.rootdist,
           crop.isgrowing,
           soil.w,
           soil.wsat,
           soil.NO3,
           soil.NH4,
           soil.layer_depth,
           soil.temp,
           ndrange=length(crop.nitrogen))
    
    KernelAbstractions.synchronize(backend)
  
end

@kernel function nuptake_crop_kernel!(PFT::PftParameters,
                                      lpjmlparams::LPJmLParams,
                                      crop_leafn::AbstractArray{T},
                                      crop_leafc::AbstractArray{T},
                                      crop_rootn::AbstractArray{T},
                                      crop_rootc::AbstractArray{T},
                                      crop_ndemand_leaf::AbstractArray{T},
                                      crop_ndemand_tot::AbstractArray{T},
                                      crop_nitrogen::AbstractArray{T},
                                      crop_vscal::AbstractArray{T},
                                      crop_rootdist::AbstractArray{T},
                                      crop_isgrowing::AbstractArray{S},
                                      soil_w::AbstractArray{M},
                                      soil_wsat::AbstractArray{M},
                                      soil_NO3::AbstractArray{M},
                                      soil_NH4::AbstractArray{M},
                                      soil_layer_depth::AbstractArray{T},
                                      soil_temp::AbstractArray{M};
                                      soil_layers = 5, # Priestley-Taylor coefficient
                                      AUTO_FERTILIZER = true
) where {T <: AbstractFloat, M <: AbstractFloat, S <: Integer}
    
    cell = @index(Global)
    
    @unpack T_0, T_m, T_r = lpjmlparams
    @unpack ncleaf, knstore, vmax_up, kNmin, KNmin = PFT

    if crop_isgrowing[cell] == 1
        NCplant = (crop_leafn[cell] + crop_rootn[cell]) / (crop_leafc[cell] + crop_rootc[cell]) # Plant's mobile nitrogen concentration
        f_NCplant = min(max(((NCplant - ncleaf.high) / (ncleaf.low - ncleaf.high)), 0), 1)

        n_uptake = zero(T)
        nsum = zero(T)

        if (crop_leafn[cell] / crop_leafc[cell]) < (ncleaf.high * (1 + knstore))
            wscaler = zero(T)
            totn = zero(T)
            NO3_up = zero(T)
            for l in 1:soil_layers
                if soil_w[l, cell] > 1.0e-7
                    wscaler = one(T)
                end
                totn = (soil_NO3[l, cell] + soil_NH4[l, cell]) * wscaler
                nuptake_temp_fcn = max((soil_temp[l, cell] - T_0) * (2 * T_m - T_0 - soil_temp[l, cell]) / (T_r - T_0) / (2 * T_m - T_0 - T_r), 0)
                if totn > 0
                    NO3_up = 2 * vmax_up * (kNmin + totn / (totn + KNmin * soil_wsat[l, cell] * soil_layer_depth[l] / 1000)) * nuptake_temp_fcn * f_NCplant * crop_rootc[cell] * crop_rootdist[l] / 1000
                end
                if NO3_up > totn
                    NO3_up = totn
                end
                n_uptake += NO3_up
                nsum += totn * crop_rootdist[l]
            end
        end
        
        if nsum == 0
            n_uptake = zero(T)
        else
            if n_uptake > (crop_ndemand_tot[cell] - crop_nitrogen[cell])
                n_uptake = crop_ndemand_tot[cell] - crop_nitrogen[cell]
            end
            if n_uptake <= 0
                n_uptake = zero(T)
            else
                crop_nitrogen[cell] += n_uptake
                wscaler = zero(T)
                for l in 1:soil_layers
                    if soil_w[l, cell] > 1.0e-7
                        wscaler = one(T)
                    end
                    soil_NO3[l, cell] -= (soil_NO3[l, cell] * wscaler * crop_rootdist[l] * n_uptake) / nsum
                    soil_NH4[l, cell] -= (soil_NH4[l, cell] * wscaler * crop_rootdist[l] * n_uptake) / nsum
                    if soil_NO3[l, cell] < 0 
                        crop_nitrogen[cell] += soil_NO3[l, cell]
                        soil_NO3[l, cell] = zero(T)
                    end
                    if soil_NH4[l, cell] < 0 
                        crop_nitrogen[cell] += soil_NH4[l, cell]
                        soil_NH4[l, cell] = zero(T)
                    end
                end
            end
        end

        ndemand_leaf_opt = crop_ndemand_leaf[cell]
        if crop_ndemand_tot[cell] > crop_nitrogen[cell]
            crop_ndemand_leaf[cell] = crop_leafn[cell]
            if ndemand_leaf_opt < 1.0e-7
                crop_vscal[cell] = one(T)
            else
                crop_vscal[cell] = min(one(T), (crop_ndemand_leaf[cell] / (ndemand_leaf_opt / (1 + knstore))))
            end
        else
            crop_vscal[cell] = one(T)
        end

    else
        crop_nitrogen[cell] = zero(T)
        crop_vscal[cell] = zero(T)
    end
end