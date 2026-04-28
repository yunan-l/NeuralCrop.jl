function ndemand_crop!(crop::Crop,
                       PFT::PftParameters,
                       photos_vmax::AbstractArray{T},
                       pet_daylength::AbstractArray{T},
                       temp::AbstractArray{T}
) where {T <: AbstractFloat}

    backend = KernelAbstractions.get_backend(crop.ndemand_tot)

    kernel = ndemand_crop_kernel!(backend)
    
    kernel(PFT,
           crop.lai, 
           crop.leafc, 
           crop.rootc, 
           crop.poolc, 
           crop.stoc, 
           crop.ndemand_leaf,
           crop.ndemand_tot,
           crop.isgrowing, 
           pet_daylength, 
           photos_vmax, 
           temp, 
           ndrange=length(crop.ndemand_tot))
    
    KernelAbstractions.synchronize(backend)
  
end

@kernel function ndemand_crop_kernel!(PFT::PftParameters,
                                      crop_lai::AbstractArray{T},
                                      crop_leafc::AbstractArray{T},
                                      crop_rootc::AbstractArray{T},
                                      crop_poolc::AbstractArray{T},
                                      crop_stoc::AbstractArray{T},
                                      crop_ndemand_leaf::AbstractArray{T},
                                      crop_ndemand_tot::AbstractArray{T},
                                      crop_isgrowing::AbstractArray{S},
                                      pet_daylength::AbstractArray{T},
                                      photos_vmax::AbstractArray{T},
                                      temp::AbstractArray{T};
                                      lpjmlparams::LPJmLParams = lpjmlparams,
                                      k_l = 0.08f0 # Priestley-Taylor coefficient
) where {T <: AbstractFloat, S <: Integer}
    
    cell = @index(Global)
    
    @unpack p, k_temp = lpjmlparams
    @unpack fpc, intc, ratio, ncleaf = PFT

    if crop_isgrowing[cell] == 1
        f_lai = zero(T)
        if crop_lai[cell] < 1.0
            f_lai = max(T(0.1), crop_lai[cell])
        else
            f_lai = exp(k_l * min(crop_lai[cell], T(7)))
        end

        # ndemand_leaf = zero(T)
        if pet_daylength[cell] == 0.0
            crop_ndemand_leaf[cell] = ncleaf.median * crop_leafc[cell]
        else
            crop_ndemand_leaf[cell] = p * T(0.02314815) / pet_daylength[cell] * photos_vmax[cell] * exp(-k_temp * (temp[cell] - T(25.0))) * f_lai + ncleaf.median * crop_leafc[cell]
        end

        nc_ratio = zero(T)
        if crop_leafc[cell] > 0.0
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