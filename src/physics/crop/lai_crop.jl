function lai_crop!(crop::Crop,
                   PFT::PftParameters
)

    backend = get_backend(crop.lai)
    
    kernel = lai_crop_kernel!(backend)
    
    kernel(crop.senescence, 
           crop.senescence0, 
           crop.lai, 
           crop.wscal, 
           crop.vscal, 
           crop.flaimax, 
           crop.laimax_adjusted, 
           crop.isgrowing, 
           PFT, 
           ndrange=length(crop.lai)
           )
    
    synchronize(backend)
  
end

@kernel function lai_crop_kernel!(crop_senescence::AbstractArray{B}, 
                                  crop_senescence0::AbstractArray{B},           
                                  crop_lai::AbstractArray{T},
                                  crop_wscal::AbstractArray{T},
                                  crop_vscal::AbstractArray{T},
                                  crop_flaimax::AbstractArray{T},
                                  crop_laimax_adjusted::AbstractArray{T},
                                  crop_isgrowing::AbstractArray{S},
                                  PFT::PftParameters
) where {T <: AbstractFloat, S <: Integer, B <: Bool}
    
    cell = @index(Global)

    @unpack sla, laimax = PFT

    if crop_isgrowing[cell] == 1
        lai0 = crop_lai[cell]
        if !crop_senescence[cell]
            crop_lai[cell] = crop_flaimax[cell] * laimax
            # scale daily LAI increment with minimum of wscal and vscal as simplest approach
            lai_inc = (crop_lai[cell] - lai0) * min(crop_wscal[cell]/T(1.5), crop_vscal[cell])
            crop_lai[cell] = lai_inc + lai0
        else
            if !crop_senescence0[cell]
                crop_laimax_adjusted[cell] = crop_lai[cell]
            end
            crop_lai[cell] = crop_flaimax[cell] * crop_laimax_adjusted[cell]
        end
        # if !crop_senescence[cell]
        #     # scale daily LAI increment with minimum of wscal and vscal as simplest approach
        #     lai_inc = (crop_lai[cell] - lai0) * min(crop_wscal[cell]/T(1.5), crop_vscal[cell])
        #     crop_lai[cell] = lai_inc + lai0
        # end
    else
        crop_lai[cell] = zero(T)
        crop_laimax_adjusted[cell] = zero(T)
    end
end


function lai_deficit!(crop::Crop,
                      PFT::PftParameters
)

    backend = get_backend(crop.lai)
    
    kernel = lai_deficit_kernel!(backend)
    
    kernel(crop.senescence, crop.biomass, crop.rootc, crop.leafc, crop.lai, crop.lai_nppdeficit, crop.isgrowing, PFT, ndrange=length(crop.lai))
    
    synchronize(backend)
  
end

@kernel function lai_deficit_kernel!(crop_senescence::AbstractArray{B},           
                                     crop_biomass::AbstractArray{T},           
                                     crop_rootc::AbstractArray{T}, 
                                     crop_leafc::AbstractArray{T},   
                                     crop_lai::AbstractArray{T},           
                                     crop_lai_nppdeficit::AbstractArray{T},
                                     crop_isgrowing::AbstractArray{S},       
                                     PFT::PftParameters
) where {T <: AbstractFloat, S <: Integer, B <: Bool}
    
    cell = @index(Global)

    @unpack sla = PFT

    if crop_isgrowing[cell] == 1
        if !crop_senescence[cell]
            if (crop_biomass[cell] - crop_rootc[cell]) >= crop_lai[cell] / sla
                crop_lai_nppdeficit[cell] = zero(T)
            else
                crop_lai_nppdeficit[cell] = crop_lai[cell] - crop_leafc[cell] * sla
                # today's lai_deficit is subtracted from tomorrow's LAI in lai_crop(), fpar_crop(), and actual_lai_crop().
                # These routines account for LAI effects on the simulation.
            end
        end
    else
        crop_lai_nppdeficit[cell] = zero(T)
    end
end