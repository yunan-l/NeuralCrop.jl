"""
crop_nitrogen!(crop, PFT, soil, photos_vmax, pet_daylength, temp)

Allocate acquired crop nitrogen among leaf, root, storage, and pool compartments.
"""
function crop_nitrogen!(crop::Crop,
                        PFT::PftParameters,
                        soil::Soil,
                        photos_vmax::AbstractArray{T},
                        pet_daylength::AbstractArray{T},
                        temp::AbstractArray{T}
) where {T <: AbstractFloat}


    ndemand_crop!(crop, PFT, photos_vmax, pet_daylength, temp)
    nuptake_crop!(crop, PFT, soil)

    launch_1D!(crop_nitrogen_kernel!,
               crop.nitrogen,
               crop.isgrowing,
               crop.leafc,
               crop.rootc,
               crop.stoc,
               crop.poolc,
               crop.leafn,
               crop.rootn,
               crop.ston,
               crop.pooln,
               PFT)

end


@kernel inbounds = true function crop_nitrogen_kernel!(
                                       crop_nitrogen::AbstractArray{T},
                                       crop_isgrowing::AbstractArray{S},
                                       crop_leafc::AbstractArray{T},
                                       crop_rootc::AbstractArray{T},
                                       crop_stoc::AbstractArray{T},
                                       crop_poolc::AbstractArray{T},
                                       crop_leafn::AbstractArray{T},
                                       crop_rootn::AbstractArray{T},
                                       crop_ston::AbstractArray{T},
                                       crop_pooln::AbstractArray{T},
                                       PFT::PftParameters
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
