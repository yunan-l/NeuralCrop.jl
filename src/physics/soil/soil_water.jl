function soil_water!(soil::Soil,
                     crop::Crop,
                     prec::AbstractArray{T}
) where {T <: AbstractFloat}

    Zygote.ignore() do
        soil.infil .= prec - crop.intercep
        infil_perc!(soil)
    end

    # compute soil water content
    soil.swc = soil.swc + soil.perc - crop.trans_layer - soil.evap

end