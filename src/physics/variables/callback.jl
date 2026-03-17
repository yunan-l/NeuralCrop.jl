# update crop variables for today
# function callback_crop!(crop::Crop,
#                         photos::Photos
# )

#     # crop.lai .*= crop.isgrowing
#     # crop.lai_nppdeficit .*= crop_cal.isgrowing
#     crop.rootc .*= crop.isgrowing
#     crop.leafc .*= crop.isgrowing
#     crop.stoc .*= crop.isgrowing
#     crop.poolc .*= crop.isgrowing
#     # crop.biomass .= crop.leafc .+ crop.rootc .+ crop.poolc .+ crop.stoc

#     # crop.npp .*= crop.isgrowing
#     crop.carbon_sum .+= (crop.leafc .+ crop.rootc .+ crop.poolc .+ crop.stoc)
#     crop.carbon_sum .*= crop.isgrowing
#     crop.biomass .+= crop.npp
#     crop.biomass .*= crop.isgrowing
#     # crop.phen .*= crop.isgrowing
#     # crop.albedo .*= crop.isgrowing
#     # crop.fpar .*= crop.isgrowing
#     # crop.apar .*= crop.isgrowing
    
#     photos.adt .*= crop.isgrowing
#     photos.adtmm .*= crop.isgrowing
#     photos.rd .*= crop.isgrowing

# end

function update_litc_tillage!(soil::Soil,
                              crop_cal::Calendar
)

    soil.litc = soil.litc .* (1 .- reshape(crop_cal.scallback, (1, :))) .* (1 .- reshape(crop_cal.hcallback, (1, :))) + 
                soil.tillage_frac * soil.litc .* reshape(crop_cal.scallback, (1, :)) +
                (soil.tillage_frac * (soil.litc .+ soil.c_input)) .* reshape(crop_cal.hcallback, (1, :)) 
end

function update_litn_tillage!(soil::Soil,
                              crop_cal::Calendar
)

    soil.litn = soil.litn .* (1 .- reshape(crop_cal.scallback, (1, :))) .* (1 .- reshape(crop_cal.hcallback, (1, :))) + 
                soil.tillage_frac * soil.litn .* reshape(crop_cal.scallback, (1, :)) +
                (soil.tillage_frac * (soil.litn .+ soil.c_input)) .* reshape(crop_cal.hcallback, (1, :)) 

end


function update_lit_winter_wheat!(soil::Soil,
                                  litch::AbstractArray{M},
                                  litnh::AbstractArray{M},
                                  crop_wtype::AbstractArray{B},
                                  hdate::AbstractArray{S},
                                  crop_cal_hcallback::AbstractArray{S},
                                  day::Int
) where {M <: AbstractFloat, B <: Bool, S <: Integer}

    hdate_callback = copy(crop_cal_hcallback)

    Zygote.ignore() do
        backend = get_backend(hdate_callback)
        kernel = update_lit_winter_wheat_kernel!(backend)
        kernel(crop_wtype,
               hdate,
               hdate_callback, 
               day,
               ndrange=length(hdate_callback)
        )
        synchronize(backend)
    end

    soil.litc = soil.litc .* (1 .- reshape(hdate_callback, (1, :))) + litch .* reshape(hdate_callback, (1, :)) 
    soil.litn = soil.litn .* (1 .- reshape(hdate_callback, (1, :))) + litnh .* reshape(hdate_callback, (1, :)) 

end


@kernel function update_lit_winter_wheat_kernel!(crop_wtype::AbstractArray{B},
                                                 hdate::AbstractArray{S},
                                                 hdate_callback::AbstractArray{S},
                                                 day::Int

) where {B <: Bool, S <: Integer}

    cell = @index(Global)

    hdate_callback[cell] = zero(S)

    if day < 365 && (crop_wtype[cell] == true) && (hdate[cell] == day)
        hdate_callback[cell] = one(S)
    end

end