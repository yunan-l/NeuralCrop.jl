"""
interception!(crop, CFT, pet_eeq, rain)

Update canopy wetness and interception evaporation for the current day.
"""
function interception!(crop,
                       CFT::CFTParameters,
                       pet_eeq::AbstractArray{T},
                       rain::AbstractArray{T};
                       lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat}

    launch_1D!(
        interception_kernel!,
        crop_fluxes(crop).water.interception,
        crop_canopy_auxiliary(crop).canopy_wet,
        crop_prognostic(crop).canopy.lai,
        crop_prognostic(crop).canopy.lai_npp_deficit,
        crop_prognostic(crop).phenology.is_growing,
        pet_eeq,
        rain,
        kernel_constant(
            crop_fluxes(crop).water.interception, (; CFT, lpjmlparams),
        ),
    )

end

@kernel inbounds = true function interception_kernel!(
                                      crop_intercep::AbstractArray{T},
                                      crop_canopy_wet::AbstractArray{T},
                                      crop_lai::AbstractArray{T},
                                      crop_lai_nppdeficit::AbstractArray{T},
                                      crop_isgrowing::AbstractArray{S},
                                      pet_eeq::AbstractArray{T},
                                      rain::AbstractArray{T},
                                      fixed_parameters,
) where {T <: AbstractFloat, S <: Integer}

    cell = @index(Global)
    parameters = kernel_value(fixed_parameters)
    CFT = parameters.CFT
    lpjmlparams = parameters.lpjmlparams

    @unpack PRIESTLEY_TAYLOR = lpjmlparams
    @unpack fpc, intc = CFT

    if crop_isgrowing[cell] == 1
        if pet_eeq[cell] < 0.0001 || fpc == 0.0
            crop_canopy_wet[cell] = zero(T)
        else
            actual_lai = max(zero(T), crop_lai[cell] - crop_lai_nppdeficit[cell])
            int_store = intc * actual_lai
            if int_store > 0.9999
                int_store = T(0.9999)
            end
            crop_canopy_wet[cell] = int_store * rain[cell] / (pet_eeq[cell] * PRIESTLEY_TAYLOR)
            if crop_canopy_wet[cell] > 0.9999
                crop_canopy_wet[cell] = T(0.9999)
            end
        end
        crop_intercep[cell] = pet_eeq[cell] * PRIESTLEY_TAYLOR * crop_canopy_wet[cell] * fpc
    else
        crop_intercep[cell] = zero(T)
        crop_canopy_wet[cell] = zero(T)
    end
end
