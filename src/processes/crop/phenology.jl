"""
phenology_crop!(crop, climbuf_V_req, CFT, temp, daylength)

Advance crop phenology, vernalization status, and harvest/senescence flags.
"""
function phenology_crop!(crop,
                         climbuf_V_req::AbstractArray{T},
                         CFT::CFTParameters,
                         temp::AbstractArray{T},
                         daylength::AbstractArray{T},
) where {T <: AbstractFloat}

    # 1D launch over cells; climbuf_V_req is used as launch reference and kernel arg #1.
    launch_1D!(
        phenology_kernel!,
        climbuf_V_req,
        crop_phenology_input(crop).phu,
        crop_prognostic(crop).phenology.vdsum,
        crop_prognostic(crop).phenology.husum,
        crop_phenology_auxiliary(crop).fphu,
        crop_canopy_auxiliary(crop).flaimax,
        crop_prognostic(crop).phenology.senescence,
        crop_prognostic(crop).phenology.senescence_previous,
        crop_prognostic(crop).phenology.harvesting,
        crop_prognostic(crop).phenology.harvesting_previous,
        crop_prognostic(crop).phenology.growing_days,
        crop_prognostic(crop).phenology.is_growing,
        crop_phenology_input(crop).winter_type,
        temp,
        daylength,
        kernel_constant(climbuf_V_req, CFT),
    )

    lai_crop!(crop, CFT)

end

@inline compute_phenology_fraction(husum::T, phu::T) where {T <: AbstractFloat} =
    phu > zero(T) ? min(one(T), husum / phu) : zero(T)

@inline compute_heat_unit_increment(temperature::T, base_temperature::T) where {T <: AbstractFloat} =
    max(zero(T), temperature - base_temperature)

@inline function compute_vernalization_increment(
    temperature::T,
    accumulated::T,
    requirement::T,
    effective,
    optimum,
) where {T <: AbstractFloat}
    accumulated < requirement || return zero(T)
    if temperature >= effective.low && temperature < optimum.low
        return (temperature - effective.low) / (optimum.low - effective.low)
    elseif temperature <= effective.high && temperature >= optimum.high
        return (effective.high - temperature) / (effective.high - optimum.high)
    elseif temperature >= optimum.low && temperature < optimum.high
        return one(T)
    end
    return zero(T)
end

@inline function compute_vernalization_factor(
    accumulated::T,
    requirement::T,
) where {T <: AbstractFloat}
    base_requirement = requirement / T(5)
    if accumulated < base_requirement
        return zero(T)
    elseif accumulated < requirement
        return max(
            zero(T), min(one(T), (accumulated - base_requirement) /
                                (requirement - base_requirement)),
        )
    end
    return one(T)
end

@inline function compute_photoperiod_factor(
    previous_fphu::T,
    senescence_fraction::T,
    sensitivity::T,
    daylength::T,
    base_daylength::T,
    saturating_daylength::T,
) where {T <: AbstractFloat}
    previous_fphu <= senescence_fraction || return one(T)
    return (one(T) - sensitivity) * min(
        one(T), max(zero(T), (daylength - base_daylength) /
                             (saturating_daylength - base_daylength)),
    ) + sensitivity
end

@inline function compute_phenology_lai_fraction(
    fphu::T,
    fphuc::T,
    flaimaxc::T,
    fphuk::T,
    flaimaxk::T,
    senescence_fraction::T,
    harvest_fraction::T,
    senescence_shape::T,
) where {T <: AbstractFloat}
    if fphu < senescence_fraction
        c = fphuc / flaimaxc - fphuc
        k = fphuk / flaimaxk - fphuk
        return fphu / (fphu + c * (c / k)^((fphuc - fphu) / (fphuk - fphuc)))
    end
    return ((one(T) - fphu) / (one(T) - senescence_fraction))^senescence_shape *
           (one(T) - harvest_fraction) + harvest_fraction
end

@kernel inbounds = true function phenology_kernel!(
                                   climbuf_V_req::AbstractArray{T},
                                   crop_phu::AbstractArray{T},
                                   crop_vdsum::AbstractArray{T},
                                   crop_husum::AbstractArray{T},
                                   crop_fphu::AbstractArray{T},
                                   crop_flaimax::AbstractArray{T},
                                   crop_senescence::AbstractArray{B},
                                   crop_senescence0::AbstractArray{B},
                                   crop_harvesting::AbstractArray{B},
                                   crop_harvesting_previous::AbstractArray{B},
                                   crop_growingdays::AbstractArray{S},
                                   crop_isgrowing::AbstractArray{S},
                                   crop_wtype::AbstractArray{B},
                                   temp::AbstractArray{T},
                                   daylength::AbstractArray{T},
                                   fixed_cft,
) where {T <: AbstractFloat, B <: Bool, S <: Integer}

    cell = @index(Global)
    CFT = kernel_value(fixed_cft)

    @unpack basetemp, tv_eff, tv_opt, fphuc, flaimaxc, fphuk, flaimaxk, fphusen, flaimaxharvest, psens, pb, ps, hlimit, sla, shapesenescencenorm = CFT

    crop_harvesting_previous[cell] = crop_harvesting[cell]
    crop_senescence0[cell] = crop_senescence[cell]

    if crop_isgrowing[cell] == 1
        crop_growingdays[cell] += 1
        if crop_husum[cell] < crop_phu[cell]
            # Daily heat-unit increment above base temperature.
            hu = compute_heat_unit_increment(temp[cell], basetemp.low)
            if crop_wtype[cell] # winter crops with vernalization requirements
                vd_inc = compute_vernalization_increment(
                    temp[cell], crop_vdsum[cell], climbuf_V_req[cell], tv_eff, tv_opt,
                )
                crop_vdsum[cell] += max(zero(T), vd_inc)
                vrf = compute_vernalization_factor(crop_vdsum[cell], climbuf_V_req[cell])
            else
                vrf = one(T)
            end

            previous_fphu = compute_phenology_fraction(crop_husum[cell], crop_phu[cell])
            prf = compute_photoperiod_factor(previous_fphu, fphusen, psens, daylength[cell], pb, ps)

            #Calculation of temperature sum (deg Cd)
            crop_husum[cell] += hu * vrf * prf

            #fraction of growing season
            crop_fphu[cell] = compute_phenology_fraction(crop_husum[cell], crop_phu[cell])

            if crop_fphu[cell] < fphusen
                crop_flaimax[cell] = compute_phenology_lai_fraction(
                    crop_fphu[cell], fphuc, flaimaxc, fphuk, flaimaxk,
                    fphusen, flaimaxharvest, shapesenescencenorm,
                )
            else
                crop_senescence[cell] = true
                crop_flaimax[cell] = compute_phenology_lai_fraction(
                    crop_fphu[cell], fphuc, flaimaxc, fphuk, flaimaxk,
                    fphusen, flaimaxharvest, shapesenescencenorm,
                )
            end

        else
            crop_harvesting[cell] = true
        end

        # `fphu` is a diagnostic cache: its process value is fully determined
        # by the prognostic heat sum and fixed PHU requirement.
        crop_fphu[cell] = compute_phenology_fraction(crop_husum[cell], crop_phu[cell])

        if(crop_growingdays[cell] == hlimit)
            crop_harvesting[cell] = true
        end
    else
        crop_vdsum[cell] = zero(T)
        crop_husum[cell] = zero(T)
        crop_fphu[cell] = zero(T)
        crop_senescence[cell] = false
        crop_harvesting[cell] = false
        crop_harvesting_previous[cell] = false
        crop_growingdays[cell] = 0
        crop_flaimax[cell] = zero(T)
    end
end
