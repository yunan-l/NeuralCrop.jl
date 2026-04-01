function crop_carbon!(photos::Photos,
                      crop::Crop,
                      PFT::PftParameters,
                      param::LPJmLParam,
                      temp::AbstractArray{T},
) where {T <: AbstractFloat} # directly translated from LPJmL

    # compute crop respiration
    respiration!(crop, PFT, param, temp, photos.agd - photos.rd)

    # compute crop carbon allocation  
    carbon_allocation!(PFT, crop, photos)
    crop.vegc = vcat(reshape(crop.rootc, (1, :)), reshape(crop.leafc, (1, :)), reshape(crop.stoc, (1, :)), reshape(crop.poolc, (1, :)))

    # crop.carbon_sum = crop.rootc .+ crop.leafc .+ crop.stoc .+ crop.poolc
    # crop.carbon_sum = crop.carbon_sum .* crop.isgrowing
end