# ifelse is more friendly to GPU parallel computing than idx and @kernel
function cultivate!(crop::Crop,
                    crop_cal::Calendar,
                    crop_sdate::AbstractArray{S},
                    param::LPJmLParam,
                    ml::Managed_land,
                    soil::Soil,
                    day::Int,
                    device
)  where {S <: Integer}

    @unpack manure_cn, nfert_split_frac, nmanure_nh4_frac= param

    Zygote.ignore() do
        # if day > 1 && day % 365 == 1
        #     crop_cal.sdate = crop_sdate[div(day, 365) + 1, :]
        # end
        # day_ = day % 365 != 0 ? day % 365 : 365
        crop.harvesting .= ifelse.(crop_cal.sdate .== day, false, crop.harvesting)
        # Update scallback and g_period
        crop_cal.scallback .= ifelse.(crop_cal.sdate .== day, 1, crop_cal.scallback)
        crop.isgrowing .= ifelse.(crop_cal.sdate .== day, 1, crop.isgrowing)
        crop_cal.scallback .= ifelse.(crop_cal.sdate .!= day, 0, crop_cal.scallback)
        fertilizer!(param, crop_cal, ml, crop, soil, day)
    end

    crop.lai = crop.lai .* (1 .- crop_cal.scallback) .+ 0.000415f0 .* crop_cal.scallback
    crop.biomass = crop.biomass .* (1 .- crop_cal.scallback) .+ 20.0f0 .* crop_cal.scallback
    init_vegc = device([8.0f0, 0.0113804f0, 0.0f0, 11.9886196f0])
    crop.vegc = crop.vegc .* (1 .- reshape(crop_cal.scallback, (1, :))) .+ init_vegc .* reshape(crop_cal.scallback, (1, :))

    # initilization of crop carbon pools and nitrogen pools
    crop.rootc = crop.rootc .* (1 .- crop_cal.scallback) .+ 8.0f0 .* crop_cal.scallback
    crop.leafc = crop.leafc .* (1 .- crop_cal.scallback) .+ 0.0113804f0 .* crop_cal.scallback
    crop.stoc = crop.stoc .* (1 .- crop_cal.scallback) .+ 0.0f0 .* crop_cal.scallback
    crop.poolc = crop.poolc .* (1 .- crop_cal.scallback) .+ 11.9886196f0 .* crop_cal.scallback

    init_nitrogen = 0.7f0 # C:N ratio of seed = 29
    crop.nitrogen = crop.nitrogen .* (1 .- crop_cal.scallback) .+ init_nitrogen * crop_cal.scallback

    # ignore it
    # soil.litc[2, :] = soil.litc[2, :] + (ml.manure * manure_cn * nfert_split_frac) * reshape(crop_cal.scallback, (1, :))
    # soil.litn[2, :] = soil.litn[2, :] + (ml.manure * (1 - nmanure_nh4_frac) * nfert_split_frac) * reshape(crop_cal.scallback, (1, :))

end