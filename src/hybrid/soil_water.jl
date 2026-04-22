function soil_water!(nn_model,
                     ps,
                     st,
                     soil::Soil,
                     crop::Crop,
                     prec::AbstractArray{T},
                     sw_n::AbstractArray{T},
                     lw_n::AbstractArray{T};
                     irrigation = false,
                     lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat}

    Zygote.ignore() do
        soil.infil .= prec - crop.intercep
        infil_perc!(soil)
    end

    #compute soil water content
    if irrigation
        # soil moisture is set to field capacity every day, without constraints on water availability
        input = vcat(mean(soil.temp, dims = 1)/30, reshape(sw_n, (1, :)), reshape(lw_n, (1, :)), mean(soil.Ks, dims = 1)/20)
        soil.swc = neural_moisture(nn_model, soil.swc, ps, st, input, soil.layer_depth, soil.perc, crop.trans_layer, soil.evap)
        # # soil.swc = neural_moisture(nn_model, soil.swc, ps, st, input, soil.layer_depth, soil.perc, crop.trans_layer)
    else
        soil.swc = soil.wfc .* soil.layer_depth
    end
end