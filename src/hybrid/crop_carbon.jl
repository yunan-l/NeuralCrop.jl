function crop_carbon!(nn_model,
                      ps,
                      st,
                      photos::Photos,
                      crop::Crop,
                      PFT::PftParameters,
                      param::LPJmLParam,
                      temp::AbstractArray{T},
                      temp_n::AbstractArray{T},
                      soil_swc::AbstractArray{M},
) where {T <: AbstractFloat, M <: AbstractFloat}

    # compute crop respiration
    Zygote.ignore() do
        respiration!(crop, PFT, param, temp, photos.agd - photos.rd)
    end

    # compute crop carbon allocation
    crop.npp = (photos.agd - photos.rd - crop.resp)
    crop.biomass = crop.biomass .+ crop.npp
    crop.biomass = crop.biomass .* crop.isgrowing
    
    # input = vcat(reshape(crop.npp/20, (1, :)), reshape(crop.lai/5, (1, :)), reshape(crop.leafn/5, (1, :)), mean(soil_swc[1:3, :], dims = 1), reshape(temp_n, (1, :))) .* reshape(crop.isgrowing, (1, :))
    input = vcat(reshape(crop.npp/20, (1, :)), reshape(crop.lai/5, (1, :)), reshape(crop.fphu, (1, :)), mean(soil_swc[1:3, :], dims = 1)) .* reshape(crop.isgrowing, (1, :))
    crop.vegc = neural_allocation(nn_model, crop.vegc, ps, st, input)
    # crop_vegc = crop.vegc ./ (sum(max.(crop.vegc, zero(T)), dims = 1) .+ T(0.00001)) .* reshape(crop.biomass, (1, :)) # for mass balance of projection
    
    crop.rootc = crop.vegc[1, :] .* crop.isgrowing
    crop.leafc = crop.vegc[2, :] .* crop.isgrowing
    crop.stoc = crop.vegc[3, :] .* crop.isgrowing
    crop.poolc = crop.vegc[4, :] .* crop.isgrowing

    # crop.carbon_sum = crop.rootc .+ crop.leafc .+ crop.stoc .+ crop.poolc
    # crop.carbon_sum = crop.carbon_sum .* crop.isgrowing
end


function crop_carbon_old!(photos::Photos,
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