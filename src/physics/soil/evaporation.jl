function evaporation!(pet_eeq::AbstractArray{T},
                      param::LPJmLParam,
                      crop::Crop,
                      soil::Soil
    
) where {T <: AbstractFloat}

    backend = KernelAbstractions.get_backend(pet_eeq)

    kernel = evaporation_kernel!(backend)
    
    kernel(pet_eeq, 
           param, 
           crop.fpar, 
           crop.trans_layer,
           crop.canopy_wet, 
           soil.swc,
           soil.wpwps,
           soil.whcs,
           soil.evap,
           soil.agtop_cover,
           soil.layer_depth,
           ndrange=length(pet_eeq))
    
    KernelAbstractions.synchronize(backend)

end

@kernel function evaporation_kernel!(pet_eeq::AbstractArray{T},
                                     param::LPJmLParam,
                                     crop_fpar::AbstractArray{T},
                                     crop_trans_layer::AbstractArray{M},
                                     crop_canopy_wet::AbstractArray{T},
                                     soil_swc::AbstractArray{M},
                                     soil_wpwps::AbstractArray{M},
                                     soil_whcs::AbstractArray{M},
                                     soil_evap::AbstractArray{M},
                                     soil_agtop_cover::AbstractArray{T},
                                     soil_layer_depth::AbstractArray{T};
                                     soil_layers = 5
) where {T <: AbstractFloat, M <: AbstractFloat}
    
    cell = @index(Global)

    @unpack PRIESTLEY_TAYLOR = param  # Priestley-Taylor coefficient
    
    soildepth_evap = param.soildepth_evap

    evap_energy = pet_eeq[cell] * PRIESTLEY_TAYLOR * max(1 - crop_fpar[cell], 0.05)
    # evap_litter = pet_eeq[cell] * PRIESTLEY_TAYLOR * (1 - crop_canopy_wet[cell]) - sum(crop_trans_layer[:, cell])

    crop_trans_layer_sum = zero(T)
    for l in 1:soil_layers
        crop_trans_layer_sum += crop_trans_layer[l, cell]
    end

    evap_ratio = zero(T)
    if evap_energy > 1.0f-5 && (pet_eeq[cell] * PRIESTLEY_TAYLOR * (1 - crop_canopy_wet[cell]) - crop_trans_layer_sum) > 1.0f-5
        # w_evap is water content in soildepth_evap that can evaporate
        w_evap = zero(T)
        whcs_evap = zero(T)

        for l in 1:soil_layers
            if soildepth_evap > 0
                fraction = min(1, soildepth_evap / soil_layer_depth[l])
                w_evap += (soil_swc[l, cell] - soil_wpwps[l, cell] - crop_trans_layer[l, cell]) * fraction
                whcs_evap += soil_whcs[l, cell] * fraction
                soildepth_evap -= soil_layer_depth[l]
            end
        end

        evap_soil = evap_energy / (1 + exp(5 - 10 * w_evap / whcs_evap)) * max(0.05, (1 - soil_agtop_cover[cell]))
        if w_evap > 0
            evap_ratio = evap_soil / w_evap
        else
            evap_ratio = zero(T)
        end
    end

    soildepth_evap = param.soildepth_evap
    for l in 1:soil_layers
        if soildepth_evap > 0
            fraction = min(1, soildepth_evap / soil_layer_depth[l])
            soil_evap[l, cell] = (soil_swc[l, cell] - soil_wpwps[l, cell]) * evap_ratio * fraction
            soildepth_evap -= soil_layer_depth[l]
        end
    end

end