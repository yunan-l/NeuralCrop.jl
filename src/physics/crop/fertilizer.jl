function fertilizer!(param::LPJmLParam,
                     crop_cal::Calendar,
                     ml::Managed_land,
                     crop::Crop,
                     soil::Soil,
                     day
)

    backend = KernelAbstractions.get_backend(crop.nfertilizer)

    kernel = fertilizer_kernel!(backend)
    
    kernel(param,
           crop_cal.sdate,
           ml.manure,
           ml.fertilizer,
           crop.nmanure,
           crop.nfertilizer,
           crop.fphu,
           soil.NO3,
           soil.NH4,
           day,
           ndrange=length(crop.nfertilizer))
    
    KernelAbstractions.synchronize(backend)

end


@kernel function fertilizer_kernel!(param::LPJmLParam,
                                    crop_cal_sdate::AbstractArray{S},
                                    ml_manure::AbstractArray{T},
                                    ml_fertilizer::AbstractArray{T},
                                    crop_nmanure::AbstractArray{T},
                                    crop_nfertilizer::AbstractArray{T},
                                    crop_fphu::AbstractArray{T},
                                    soil_NO3::AbstractArray{M},
                                    soil_NH4::AbstractArray{M},
                                    day::Integer;
                                    nitrogen_is_unlimited = true
) where {T <: AbstractFloat, M <: AbstractFloat, S <: Integer}
    
    cell = @index(Global)

    @unpack nmanure_nh4_frac, nfert_split_frac, nfert_no3_frac = param

    if crop_cal_sdate[cell] == day
        soil_NH4[1, cell] += ml_manure[cell] * nmanure_nh4_frac * nfert_split_frac
        crop_nmanure[cell] = ml_manure[cell] * (1 - nfert_split_frac)

        soil_NO3[1, cell] += ml_fertilizer[cell] * nfert_no3_frac * nfert_split_frac
        soil_NH4[1, cell] += ml_fertilizer[cell] * (1 - nfert_no3_frac) * nfert_split_frac
        crop_nfertilizer[cell] = ml_fertilizer[cell] * (1 - nfert_split_frac)
    end

    if crop_fphu[cell] > 0.25 && crop_nfertilizer[cell] > 0
        soil_NO3[1, cell] += crop_nfertilizer[cell]  * nfert_no3_frac
        soil_NH4[1, cell] += crop_nfertilizer[cell] * (1 - nfert_no3_frac)
        crop_nfertilizer[cell] = zero(T)
    end

    if crop_fphu[cell] > 0.25 && crop_nmanure[cell] > 0
        soil_NH4[1, cell] += crop_nmanure[cell] * nmanure_nh4_frac
        crop_nmanure[cell] = zero(T)
    end

    if nitrogen_is_unlimited
        soil_NO3[1, cell] += T(1000)
        soil_NH4[1, cell] += T(1000)
    end

end