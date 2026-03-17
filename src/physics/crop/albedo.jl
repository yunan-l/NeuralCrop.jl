
function albedo!(PFT::PftParameters,
                 crop::Crop,
                 pet_albedo::AbstractArray{T},
                 soil_albedo = 0.3f0  # Albedo of bare soil (0-1). Should be soil and soil moisture dependent */
) where {T <: AbstractFloat}

    crop_albedo!(PFT, crop)

    pet_albedo .= crop.albedo .+ (1 .- PFT.fpc * crop.isgrowing) * soil_albedo

end


function crop_albedo!(PFT::PftParameters,
                      crop::Crop
)
    @unpack albedo_leaf, albedo_litter, fpc = PFT
    
    albedo_green_leaves = fpc * crop.phen * albedo_leaf
    
    # albedo of PFT without green foliage (litter background albedo)
  
    albedo_brown_litter = fpc * (1 .- crop.phen) * albedo_litter
    
    crop.albedo .= albedo_green_leaves .+ albedo_brown_litter
    
end
