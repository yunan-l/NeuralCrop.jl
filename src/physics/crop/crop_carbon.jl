function crop_carbon!(photos::Photos,
                      crop::Crop,
                      PFT::PftParameters,
                      temp::AbstractArray{T};
                      lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat} # directly translated from LPJmL

    # compute crop respiration
    respiration!(crop, PFT, temp, photos.agd - photos.rd)

    # compute crop carbon allocation  
    carbon_allocation!(PFT, crop, photos)
    crop.vegc = vcat(reshape(crop.rootc, (1, :)), reshape(crop.leafc, (1, :)), reshape(crop.stoc, (1, :)), reshape(crop.poolc, (1, :)))

end