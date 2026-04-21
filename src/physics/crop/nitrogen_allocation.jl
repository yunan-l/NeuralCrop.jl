function crop_nitrogen!(crop::Crop,
                        PFT::PftParameters,
                        soil::Soil,
                        photos_vmax::AbstractArray{T},
                        pet_daylength::AbstractArray{T},
                        temp::AbstractArray{T}
) where {T <: AbstractFloat}


    ndemand_crop!(crop, PFT, photos_vmax, pet_daylength, temp)
    nuptake_crop!(crop, PFT, soil)

    backend = KernelAbstractions.get_backend(crop.nitrogen)

    kernel = crop_nitrogen_kernel!(backend)
    
    kernel(PFT,
           crop.isgrowing,
           crop.nitrogen,
           crop.leafc,
           crop.rootc,
           crop.stoc,
           crop.poolc,
           crop.leafn,
           crop.rootn,
           crop.ston,
           crop.pooln,
           ndrange=length(crop.nitrogen))
    
    KernelAbstractions.synchronize(backend)

end


@kernel function crop_nitrogen_kernel!(PFT::PftParameters,
                                       crop_isgrowing::AbstractArray{S},
                                       crop_nitrogen::AbstractArray{T},
                                       crop_leafc::AbstractArray{T},
                                       crop_rootc::AbstractArray{T},
                                       crop_stoc::AbstractArray{T},
                                       crop_poolc::AbstractArray{T},
                                       crop_leafn::AbstractArray{T},
                                       crop_rootn::AbstractArray{T},
                                       crop_ston::AbstractArray{T},
                                       crop_pooln::AbstractArray{T}
) where {T <: AbstractFloat, S <: Integer}

     cell = @index(Global)

     @unpack ratio = PFT
 
     if (crop_isgrowing[cell] == 1) && (crop_nitrogen[cell] > 0)
          nominator = crop_nitrogen[cell] * (crop_leafc[cell] * ratio.root * ratio.sto * ratio.pool +
                      crop_rootc[cell] * ratio.sto * ratio.pool +
                      crop_stoc[cell] * ratio.root * ratio.pool +
                      crop_poolc[cell] * ratio.root * ratio.sto)

          a = (crop_leafc[cell] * crop_nitrogen[cell] * ratio.root * ratio.sto * ratio.pool +
               crop_leafc[cell] * crop_rootn[cell] * ratio.root * ratio.sto * ratio.pool +
               crop_leafc[cell] * crop_ston[cell] * ratio.root * ratio.sto * ratio.pool +
               crop_leafc[cell] * crop_pooln[cell] * ratio.root * ratio.sto * ratio.pool -
               crop_rootc[cell] * crop_leafn[cell] * ratio.sto * ratio.pool -
               crop_stoc[cell] * crop_leafn[cell] * ratio.root * ratio.pool -
               crop_poolc[cell] * crop_leafn[cell] * ratio.root * ratio.sto) / nominator

          b = (-crop_leafc[cell] * crop_rootn[cell] * ratio.root * ratio.sto * ratio.pool +
               crop_rootc[cell] * crop_nitrogen[cell] * ratio.sto * ratio.pool +
               crop_rootc[cell] * crop_leafn[cell] * ratio.sto * ratio.pool +
               crop_rootc[cell] * crop_ston[cell] * ratio.sto * ratio.pool +
               crop_rootc[cell] * crop_pooln[cell] * ratio.sto * ratio.pool -
               crop_stoc[cell] * crop_rootn[cell] * ratio.root * ratio.pool -
               crop_poolc[cell] * crop_rootn[cell] * ratio.root * ratio.sto) / nominator

          c = (-crop_leafc[cell] * crop_ston[cell] * ratio.root * ratio.sto * ratio.pool -
               crop_rootc[cell] * crop_ston[cell] * ratio.sto * ratio.pool +
               crop_stoc[cell] * crop_nitrogen[cell] * ratio.root * ratio.pool +
               crop_stoc[cell] * crop_leafn[cell] * ratio.root * ratio.pool +
               crop_stoc[cell] * crop_rootn[cell] * ratio.root * ratio.pool +
               crop_stoc[cell] * crop_pooln[cell] * ratio.root * ratio.pool -
               crop_poolc[cell] * crop_ston[cell] * ratio.root * ratio.sto) / nominator

          d = (-crop_leafc[cell] * crop_pooln[cell] * ratio.root * ratio.sto * ratio.pool -
               crop_rootc[cell] * crop_pooln[cell] * ratio.sto * ratio.pool -
               crop_stoc[cell] * crop_pooln[cell] * ratio.root * ratio.pool +
               crop_poolc[cell] * crop_nitrogen[cell] * ratio.root * ratio.sto +
               crop_poolc[cell] * crop_leafn[cell] * ratio.root * ratio.sto +
               crop_poolc[cell] * crop_rootn[cell] * ratio.root * ratio.sto +
               crop_poolc[cell] * crop_ston[cell] * ratio.root * ratio.sto) / nominator

          crop_leafn[cell] += a * crop_nitrogen[cell]
          crop_rootn[cell] += b * crop_nitrogen[cell]
          crop_ston[cell] += c * crop_nitrogen[cell]
          crop_pooln[cell] += d * crop_nitrogen[cell]
     else
          crop_leafn[cell] = zero(T)
          crop_rootn[cell] = zero(T)
          crop_ston[cell] = zero(T)
          crop_pooln[cell] = zero(T)

     end
end

function crop_nitrogen_old!(crop::Crop,
                            PFT::PftParameters,
                            soil::Soil,
                            param::LPJmLParams,
                            photos_vmax::AbstractArray{T},
                            pet_daylength::AbstractArray{T},
                            temp::AbstractArray{T}
) where {T <: AbstractFloat}


    ndemand_crop!(crop, PFT, photos_vmax, pet_daylength, temp)
    nuptake_crop!(crop, PFT, soil)

    backend = KernelAbstractions.get_backend(crop.nitrogen)

    kernel = crop_nitrogen_old_kernel!(backend)
    
    kernel(PFT,
           crop.isgrowing,
           crop.nitrogen,
           crop.leafc,
           crop.rootc,
           crop.stoc,
           crop.poolc,
           crop.leafn,
           crop.rootn,
           crop.ston,
           crop.pooln,
           ndrange=length(crop.nitrogen))
    
    KernelAbstractions.synchronize(backend)

end


@kernel function crop_nitrogen_old_kernel!(PFT::PftParameters,
                                           crop_isgrowing::AbstractArray{S},
                                           crop_nitrogen::AbstractArray{T},
                                           crop_leafc::AbstractArray{T},
                                           crop_rootc::AbstractArray{T},
                                           crop_stoc::AbstractArray{T},
                                           crop_poolc::AbstractArray{T},
                                           crop_leafn::AbstractArray{T},
                                           crop_rootn::AbstractArray{T},
                                           crop_ston::AbstractArray{T},
                                           crop_pooln::AbstractArray{T}
) where {T <: AbstractFloat, S <: Integer}

     cell = @index(Global)

     @unpack ratio = PFT
 
     if (crop_isgrowing[cell] == 1) && (crop_nitrogen[cell] > 0)
          nominator = crop_nitrogen[cell] * (crop_leafc[cell] * ratio.root * ratio.sto * ratio.pool +
                      crop_rootc[cell] * ratio.sto * ratio.pool +
                      crop_stoc[cell] * ratio.root * ratio.pool +
                      crop_poolc[cell] * ratio.root * ratio.sto)

          a = (crop_leafc[cell] * crop_nitrogen[cell] * ratio.root * ratio.sto * ratio.pool +
               crop_leafc[cell] * crop_rootn[cell] * ratio.root * ratio.sto * ratio.pool +
               crop_leafc[cell] * crop_ston[cell] * ratio.root * ratio.sto * ratio.pool +
               crop_leafc[cell] * crop_pooln[cell] * ratio.root * ratio.sto * ratio.pool -
               crop_rootc[cell] * crop_leafn[cell] * ratio.sto * ratio.pool -
               crop_stoc[cell] * crop_leafn[cell] * ratio.root * ratio.pool -
               crop_poolc[cell] * crop_leafn[cell] * ratio.root * ratio.sto) / nominator

          b = (-crop_leafc[cell] * crop_rootn[cell] * ratio.root * ratio.sto * ratio.pool +
               crop_rootc[cell] * crop_nitrogen[cell] * ratio.sto * ratio.pool +
               crop_rootc[cell] * crop_leafn[cell] * ratio.sto * ratio.pool +
               crop_rootc[cell] * crop_ston[cell] * ratio.sto * ratio.pool +
               crop_rootc[cell] * crop_pooln[cell] * ratio.sto * ratio.pool -
               crop_stoc[cell] * crop_rootn[cell] * ratio.root * ratio.pool -
               crop_poolc[cell] * crop_rootn[cell] * ratio.root * ratio.sto) / nominator

          c = (-crop_leafc[cell] * crop_ston[cell] * ratio.root * ratio.sto * ratio.pool -
               crop_rootc[cell] * crop_ston[cell] * ratio.sto * ratio.pool +
               crop_stoc[cell] * crop_nitrogen[cell] * ratio.root * ratio.pool +
               crop_stoc[cell] * crop_leafn[cell] * ratio.root * ratio.pool +
               crop_stoc[cell] * crop_rootn[cell] * ratio.root * ratio.pool +
               crop_stoc[cell] * crop_pooln[cell] * ratio.root * ratio.pool -
               crop_poolc[cell] * crop_ston[cell] * ratio.root * ratio.sto) / nominator

          d = (-crop_leafc[cell] * crop_pooln[cell] * ratio.root * ratio.sto * ratio.pool -
               crop_rootc[cell] * crop_pooln[cell] * ratio.sto * ratio.pool -
               crop_stoc[cell] * crop_pooln[cell] * ratio.root * ratio.pool +
               crop_poolc[cell] * crop_nitrogen[cell] * ratio.root * ratio.sto +
               crop_poolc[cell] * crop_leafn[cell] * ratio.root * ratio.sto +
               crop_poolc[cell] * crop_rootn[cell] * ratio.root * ratio.sto +
               crop_poolc[cell] * crop_ston[cell] * ratio.root * ratio.sto) / nominator

          crop_leafn[cell] = a * crop_nitrogen[cell]
          crop_rootn[cell] = b * crop_nitrogen[cell]
          crop_ston[cell] = c * crop_nitrogen[cell]
          crop_pooln[cell] = d * crop_nitrogen[cell]
     else
          crop_leafn[cell] = zero(T)
          crop_rootn[cell] = zero(T)
          crop_ston[cell] = zero(T)
          crop_pooln[cell] = zero(T)

     end
end