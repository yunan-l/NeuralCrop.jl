"""
lai_crop!(crop, CFT)

Update leaf-area index from phenology and carbon state.
"""
function lai_crop!(crop,
                   CFT::CFTParameters
)

    launch_1D!(
        lai_crop_kernel!,
        crop_prognostic(crop).canopy.lai,
        crop_prognostic(crop).canopy.lai_previous_potential,
        crop_prognostic(crop).phenology.senescence,
        crop_prognostic(crop).phenology.senescence_previous,
        crop_prognostic(crop).water.sufficiency,
        crop_prognostic(crop).nitrogen.sufficiency,
        crop_canopy_auxiliary(crop).flaimax,
        crop_prognostic(crop).canopy.laimax_adjusted,
        crop_prognostic(crop).phenology.is_growing,
        kernel_constant(crop_prognostic(crop).canopy.lai, CFT),
    )

end

@kernel inbounds = true function lai_crop_kernel!(
                                  crop_lai::AbstractArray{T},
                                  crop_lai_previous_potential::AbstractArray{T},
                                  crop_senescence::AbstractArray{B},
                                  crop_senescence0::AbstractArray{B},
                                  crop_wscal::AbstractArray{T},
                                  crop_vscal::AbstractArray{T},
                                  crop_flaimax::AbstractArray{T},
                                  crop_laimax_adjusted::AbstractArray{T},
                                  crop_isgrowing::AbstractArray{S},
                                  fixed_cft,
) where {T <: AbstractFloat, S <: Integer, B <: Bool}

    cell = @index(Global)
    CFT = kernel_value(fixed_cft)

    @unpack sla, laimax = CFT

    if crop_isgrowing[cell] == 1
        lai0 = crop_lai[cell]
        if !crop_senescence[cell]
            potential_lai = crop_flaimax[cell] * laimax
            # LPJmL's `lai000` stores the previous *potential* LAI, distinct
            # from the actual leaf area retained after water/N limitation.
            # Keeping that state prevents a newly sown winter crop from losing
            # its seed LAI while vernalization keeps potential LAI at zero.
            lai_inc = (potential_lai - crop_lai_previous_potential[cell]) *
                min(crop_wscal[cell], crop_vscal[cell])
            crop_lai_previous_potential[cell] = potential_lai
            crop_lai[cell] = lai_inc + lai0
        else
            if !crop_senescence0[cell]
                crop_laimax_adjusted[cell] = crop_lai[cell]
            end
            crop_lai[cell] = crop_flaimax[cell] * crop_laimax_adjusted[cell]
        end
    else
        crop_lai[cell] = zero(T)
        crop_lai_previous_potential[cell] = zero(T)
        crop_laimax_adjusted[cell] = zero(T)
    end
end
