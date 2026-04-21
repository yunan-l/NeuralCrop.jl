function hybrid_photos_C3!(nn_model,
                           ps,
                           st,
                           PFT::PftParameters,
                           photos::Photos,
                           crop::Crop,
                           pet_daylength::AbstractArray{T},
                           soil_swc::AbstractArray{M},
                           temp_n::AbstractArray{T},
                           temp::AbstractArray{T},
                           co2
) where {T <: AbstractFloat, M <: AbstractFloat}

    # compute photosynthesis variables: lambda and vcmax
    input = vcat(reshape(pet_daylength/24, (1, :)), mean(soil_swc[1:3, :], dims = 1), reshape(temp_n, (1, :))) .* reshape(crop.isgrowing, (1, :))
    photos.lambda = neural_lambda(nn_model.lambda, ps.ps_lambda, st.st_lambda, input)

    input = vcat(reshape(pet_daylength/24, (1, :)), reshape(crop.apar/10^7, (1, :)), reshape(crop.leafn/5, (1, :)), reshape(temp_n, (1, :))) .* reshape(crop.isgrowing, (1, :))
    photos.vmax = neural_vmax(nn_model.vmax, ps.ps_vmax, st.st_vmax, input)

    # compute photosynthesis
    photosynthesis_C3!(PFT, photos, crop.apar, pet_daylength, temp, co2)

    # photos.adt = photos.adt .* crop.isgrowing
    # photos.adtmm = photos.adtmm .* crop.isgrowing
    # photos.rd = photos.rd .* crop.isgrowing
end


function hybrid_photos_C4!(nn_model,
                           ps,
                           st,
                           PFT::PftParameters,
                           photos::Photos,
                           crop::Crop,
                           pet_daylength::AbstractArray{T},
                           soil_swc::AbstractArray{M},
                           temp_n::AbstractArray{T},
                           temp::AbstractArray{T}
) where {T <: AbstractFloat, M <: AbstractFloat}


    # compute photosynthesis variables: lambda and vcmax
    input = vcat(reshape(pet_daylength/24, (1, :)), mean(soil_swc[1:3, :], dims = 1), reshape(temp_n, (1, :))) .* reshape(crop.isgrowing, (1, :))
    photos.lambda = neural_lambda(nn_model.lambda, ps.ps_lambda, st.st_lambda, input)

    input = vcat(reshape(pet_daylength/24, (1, :)), reshape(crop.apar/10^7, (1, :)), reshape(crop.leafn/5, (1, :)), reshape(temp_n, (1, :))) .* reshape(crop.isgrowing, (1, :))
    photos.vmax = neural_vmax(nn_model.vmax, ps.ps_vmax, st.st_vmax, input)

    # compute photosynthesis
    photosynthesis_C4!(PFT, photos, crop.apar, pet_daylength, temp)

    # photos.adt = photos.adt .* crop.isgrowing
    # photos.adtmm = photos.adtmm .* crop.isgrowing
    # photos.rd = photos.rd .* crop.isgrowing
end
