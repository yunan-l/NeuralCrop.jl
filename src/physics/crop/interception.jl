function interception!(crop::Crop,
                       param::LPJmLParam,
                       PFT::PftParameters,
                       pet_eeq::AbstractArray{T},
                       rain::AbstractArray{T}
) where {T <: AbstractFloat}

    backend = KernelAbstractions.get_backend(crop.intercep)

    kernel = interception_kernel!(backend)
    
    kernel(param, PFT, crop.canopy_wet, crop.lai, crop.intercep, crop.isgrowing, pet_eeq, rain, ndrange=length(crop.intercep))
    
    KernelAbstractions.synchronize(backend)
  
end

@kernel function interception_kernel!(param::LPJmLParam,
                                      PFT::PftParameters,
                                      crop_canopy_wet::AbstractArray{T},
                                      crop_lai::AbstractArray{T},
                                      crop_intercep::AbstractArray{T},
                                      crop_isgrowing::AbstractArray{S},
                                      pet_eeq::AbstractArray{T},
                                      rain::AbstractArray{T}
) where {T <: AbstractFloat, S <: Integer}
    
    cell = @index(Global)

    @unpack PRIESTLEY_TAYLOR = param
    @unpack fpc, intc = PFT

    if crop_isgrowing[cell] == 1
        if pet_eeq[cell] < 0.0001 || fpc == 0.0
            crop_canopy_wet[cell] = zero(T)
        else
            int_store = intc * crop_lai[cell]
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
    end
end