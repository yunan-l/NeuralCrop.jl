function soil_water!(nn_model,
                     ps,
                     st,
                     parm::LPJmLParam,
                     soil::Soil,
                     crop::Crop,
                     prec::AbstractArray{T},
                     sw_n::AbstractArray{T},
                     lw_n::AbstractArray{T}
) where {T <: AbstractFloat}

    Zygote.ignore() do
        soil.infil .= prec - crop.intercep
        infil_perc!(parm, soil)
    end

    #compute soil water content
    input = vcat(mean(soil.temp, dims = 1)/30, reshape(sw_n, (1, :)), reshape(lw_n, (1, :)), mean(soil.Ks, dims = 1)/20)
    soil.swc = neural_moisture(nn_model, soil.swc, ps, st, input, soil.layer_depth, soil.perc, crop.trans_layer, soil.evap)
    # # soil.swc = neural_moisture(nn_model, soil.swc, ps, st, input, soil.layer_depth, soil.perc, crop.trans_layer)

end