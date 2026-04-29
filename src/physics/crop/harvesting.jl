function harvest_crop!(crop_cal::Calendar,
                       crop::Crop,
                       soil::Soil,
                       output::Output,
                       residue_frac::AbstractArray{T},
                       device,
                       cell_size,
                       day::Int
) where {T <: AbstractFloat}

    Zygote.ignore() do
        # update hcallback and g_period
        crop_cal.hdate .= ifelse.((crop.harvesting0 .== false) .& (crop.harvesting .== true), day, crop_cal.hdate)
        crop_cal.hcallback .= ifelse.((crop.harvesting0 .== false) .& (crop.harvesting .== true), 1, 0)
        crop.isgrowing .= ifelse.((crop.harvesting0 .== false) .& (crop.harvesting .== true), 0, crop.isgrowing)
        # Update crop variables
        crop.yield .= ifelse.(((crop.harvesting0 .== false) .& (crop.harvesting .== true)), crop.stoc, crop.yield)

        soil.c_input .= vcat(reshape((crop.leafc .+ crop.poolc) .* residue_frac, (1, :)), device(zeros(Float32, (1, cell_size))), reshape(crop.rootc, (1, :))) .* reshape(crop_cal.hcallback, (1, :))
        soil.n_input .= vcat(reshape((crop.leafn .+ crop.pooln) .* residue_frac, (1, :)), device(zeros(Float32, (1, cell_size))), reshape(crop.rootn, (1, :))) .* reshape(crop_cal.hcallback, (1, :))
        # idx = ((crop.harvesting0 .== true) .& (crop.harvesting .== true)) .| ((crop.harvesting0 .== true) .& (crop.harvesting .== false)) .| ((crop.harvesting0 .== false) .& (crop.harvesting .== false))
        # crop_cal.hcallback[idx] .= 0
        # crop_cal.hcallback .= ifelse.(((crop.harvesting0 .== true) .& (crop.harvesting .== true)) .| ((crop.harvesting0 .== true) .& (crop.harvesting .== false)) .| ((crop.harvesting0 .== false) .& (crop.harvesting .== false)), 0, crop_cal.hcallback)
    end

    # update harvesting variables
    output.growing_mask = vcat(output.growing_mask, reshape(crop.isgrowing, (1, :)))
    output.harvesting_mask = vcat(output.harvesting_mask, reshape(crop_cal.hcallback, (1, :)))
    output.stoc = vcat(output.stoc, reshape(crop.stoc, (1, :)))
    output.fphu = vcat(output.fphu, reshape(crop.fphu, (1, :)))
    if day == 365
        Zygote.ignore() do
            crop_cal.harvesting_year .= ifelse.(crop.yield .!= 0.0f0, 1, 0)
            output.yield = vcat(output.yield, reshape(max.(crop.yield, 0.0f0), (1, :)))
            crop.yield .= 0.0f0
        end
        output.harvesting_year = vcat(output.harvesting_year, reshape(crop_cal.harvesting_year, (1, :)))
    end
    crop.vegc = crop.vegc .* (1 .- reshape(crop_cal.hcallback, (1, :)))
    crop.rootc = crop.rootc .* (1 .- crop_cal.hcallback)
    crop.leafc = crop.leafc .* (1 .- crop_cal.hcallback)
    crop.stoc = crop.stoc .* (1 .- crop_cal.hcallback)
    crop.poolc = crop.poolc .* (1 .- crop_cal.hcallback)

end