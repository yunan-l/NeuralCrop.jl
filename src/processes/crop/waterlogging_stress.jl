
"""
waterlogging_stress!(PFT, crop, pet)

Calculate waterlogging stress.
"""
# still under development, not tested in the current version

function waterlogging_stress!(crop::Crop,
                              soil::Soil,
                              photos_agd::AbstractArray{T}
) where {T <: AbstractFloat}
 
    # kernel function parameters
    kernel_params = (rootmoist_layers = 3, waterlogging_threshold = 0.25f0, waterlogging_recovery = 0.25f0,  waterlogging_tolerance = 3.0f0, waterlogging_maxdamage = 0.9f0, eps = 1f-7)
    
    launch_1D!(
        waterlogging_stress_kernel!,
        crop.waterlogging_days,
        crop.waterlogging_stress,
        crop.rootdist,
        crop.isgrowing,
        soil.wsats,
        soil.wpwps,
        soil.whcs,
        soil.w,
        soil.w_fw,
        kernel_params
    )

    photos_agd .*= crop.waterlogging_stress
    crop.wdf .*= crop.waterlogging_stress
end
    
    
@kernel inbounds = true function waterlogging_stress_kernel!(
                                       crop_waterlogging_days::AbstractArray{T},
                                       crop_waterlogging_stress::AbstractArray{T},
                                       crop_rootdist::AbstractArray{T},
                                       crop_isgrowing::AbstractArray{S},
                                       soil_wsats::AbstractArray{M},
                                       soil_wpwps::AbstractArray{M},
                                       soil_whcs::AbstractArray{M},
                                       soil_w::AbstractArray{M},
                                       soil_w_fw::AbstractArray{M},
                                       kernel_params 
) where {T <: AbstractFloat, M <: AbstractFloat, S <: Integer}
    
    cell = @index(Global)

    @unpack rootmoist_layers, waterlogging_threshold, waterlogging_recovery, waterlogging_tolerance, waterlogging_maxdamage, eps = kernel_params
    
    if crop_isgrowing[cell] == 1
        exposure = zero(T)
        root_sum = zero(T)

        for l in 1:rootmoist_layers
            drainable_water = soil_wsats[l, cell] - soil_wpwps[l, cell] - soil_whcs[l, cell]
            
            # Fraction of pore space between field capacity and saturation that is
            # occupied by liquid water. This is comparable across soil textures.
            if drainable_water > eps && crop_rootdist[l] > 0
                excess_saturation = (soil_w[l, cell] * soil_whcs[l, cell] + soil_w_fw[l, cell] - soil_whcs[l, cell]) / drainable_water
                layer_exposure = min(max((excess_saturation - waterlogging_threshold) / (one(T) - waterlogging_threshold), zero(T)), one(T))
                exposure += crop_rootdist[l] * layer_exposure
                root_sum += crop_rootdist[l]
            end
        end

        if root_sum > eps
            exposure = exposure / root_sum
        end

        crop_waterlogging_days[cell] = max(crop_waterlogging_days[cell] + exposure - waterlogging_recovery * (one(T) - exposure), zero(T))
        severity = min(crop_waterlogging_days[cell] / waterlogging_tolerance, one(T))
        
        crop_waterlogging_stress[cell] = max(one(T) - waterlogging_maxdamage * severity, zero(T))
    else
        crop_waterlogging_days[cell] = zero(T)
        crop_waterlogging_stress[cell] = one(T)
    end
end