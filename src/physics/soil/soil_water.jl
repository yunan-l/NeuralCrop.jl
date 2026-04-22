function soil_water!(soil::Soil,
                     crop::Crop,
                     prec::AbstractArray{T};
                     irrigation = false
) where {T <: AbstractFloat}

    Zygote.ignore() do
        soil.infil .= prec - crop.intercep
        infil_perc!(soil)
    end

    # compute soil water content
    if irrigation
        # soil moisture is set to field capacity every day, without constraints on water availability
        soil.swc = soil.wfc .* soil.layer_depth 
    else
        soil.swc = soil.swc + soil.perc - crop.trans_layer - soil.evap
    end
end