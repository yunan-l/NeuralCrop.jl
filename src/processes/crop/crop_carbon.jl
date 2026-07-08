"""
crop_carbon!(photos, PFT, crop, pet, soil, temp, co2)

Run the daily crop carbon process chain: respiration, allocation, and phenology coupling.
"""
function crop_carbon!(photos::Photos,
                      crop::Crop,
                      output::Output,
                      PFT::PftParameters,
                      temp::AbstractArray{T}
) where {T <: AbstractFloat} # directly translated from LPJmL

    # compute crop respiration
    respiration!(crop, PFT, temp, photos.agd - photos.rd)

    # compute crop carbon allocation  
    carbon_allocation!(PFT, crop, photos)
    crop.vegc = vcat(reshape(crop.rootc, (1, :)), reshape(crop.leafc, (1, :)), reshape(crop.stoc, (1, :)), reshape(crop.poolc, (1, :)))

    output.npp = vcat(output.npp, reshape(crop.npp, (1, :)))
    output.lai = vcat(output.lai, reshape(crop.lai, (1, :)))
    output.fphu = vcat(output.fphu, reshape(crop.fphu, (1, :)))
    output.biomass = vcat(output.biomass, reshape(crop.biomass, (1, :)))
          
end
