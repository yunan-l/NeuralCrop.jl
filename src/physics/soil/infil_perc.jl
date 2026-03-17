function infil_perc!(parm::LPJmLParam,
                     soil::Soil  
)

    backend = get_backend(soil.infil)

    kernel = infil_perc_kernel!(backend)

    kernel(soil.infil,           
           soil.w,       
           soil.whcs,       
           soil.w_fw,
           soil.wsat,
           soil.wsats,
           soil.wpwp,
           soil.wpwps,               
           soil.w_influx,
           soil.w_outflux,
           soil.Ks,
           soil.srunoff,
           soil.lrunoff,
           soil.outflux_f,
           soil.perc,
           soil.agtop_cover,
           soil.beta_soil,
           soil.NO3,
           soil.layer_depth,
           parm, 
           ndrange=length(soil.infil)
    )

    synchronize(backend)

end

@kernel function infil_perc_kernel!(infil::AbstractArray{T},           
                                    soil_w::AbstractArray{M},       
                                    soil_whcs::AbstractArray{M},       
                                    soil_w_fw::AbstractArray{M},
                                    soil_wsat::AbstractArray{M},        
                                    soil_wsats::AbstractArray{M},
                                    soil_wpwp::AbstractArray{M},    
                                    soil_wpwps::AbstractArray{M},               
                                    soil_w_influx::AbstractArray{M},
                                    soil_w_outflux::AbstractArray{M},
                                    soil_Ks::AbstractArray{M},
                                    soil_srunoff::AbstractArray{T},
                                    soil_lrunoff::AbstractArray{M},
                                    soil_outflux_f::AbstractArray{T},
                                    soil_perc::AbstractArray{M},
                                    soil_agtop_cover::AbstractArray{T},
                                    soil_beta_soil::AbstractArray{M},
                                    soil_NO3::AbstractArray{M},
                                    soil_layer_depth::AbstractArray{T},
                                    parm::LPJmLParam;
                                    soil_layers = 5,
                                    anion_excl = 0.3,
                                    NPERCO = 0.4
) where {T <: AbstractFloat, M <: AbstractFloat}
    
    cell = @index(Global)

    @unpack soil_infil, soil_infil_litter, percthres = parm

    freewater = zero(T)
    soil_srunoff[cell] = zero(T)
    soil_outflux_f[cell] = zero(T)
    for l in 1:soil_layers
        freewater += soil_w_fw[l, cell]
        if (soil_w[l, cell] / soil_whcs[l, cell]) > 1
            freewater += (soil_w[l, cell] - 1) * soil_whcs[l, cell]
        end
        soil_lrunoff[l, cell] = zero(T)
        soil_w_influx[l, cell] = zero(T)
        soil_w_outflux[l, cell] = zero(T)
    end

    soil_infil *= (1 + soil_agtop_cover[cell] * soil_infil_litter)
    
    influx = zero(T)
    
    iter = 0
    while (infil[cell] > 1.0f-5 || freewater > 1.0f-5) && iter < 500 ## avoid infinite loop
        iter += 1
        NO3perc_ly = zero(T)
        freewater = zero(T)
        slug = min(4, infil[cell])
        infil[cell] -= slug
        
        # Calculate influx to first soil layer
        if 1 - (soil_w[1, cell] * soil_whcs[1, cell] + soil_w_fw[1, cell]) / (soil_wsats[1, cell] - soil_wpwps[1, cell]) >= 0
            influx = slug * ((1 - (soil_w[1, cell] * soil_whcs[1, cell] + soil_w_fw[1, cell]) / (soil_wsats[1, cell] - soil_wpwps[1, cell])) ^ (1 / soil_infil))
            soil_w_influx[1, cell] += influx
        else
            influx = zero(T)
            soil_w_influx[1, cell] += influx
        end
        srunoff = slug-influx
        soil_srunoff[cell] += slug - influx # surface runoff used for leaching

        for l in 1:soil_layers
            lrunoff = zero(T)
            soil_w[l, cell] += (soil_w_fw[l, cell] + influx) / soil_whcs[l, cell]
            soil_w_fw[l, cell] = zero(T)

            # Handle lateral runoff of water above saturation
            if (soil_w[l, cell] * soil_whcs[l, cell]) > (soil_layer_depth[l] * (soil_wsat[l, cell] - soil_wpwp[l, cell]))
                grunoff = (soil_w[l, cell] * soil_whcs[l, cell]) - (soil_layer_depth[l] * (soil_wsat[l, cell] - soil_wpwp[l, cell]))
                soil_w[l, cell] -= grunoff / soil_whcs[l, cell]
                soil_lrunoff[l, cell] += grunoff
                lrunoff += grunoff
            end
            
            # Additional saturation check
            if (soil_wpwps[l, cell] + soil_w[l, cell] * soil_whcs[l, cell]) > soil_wsats[l, cell]
                grunoff = (soil_wpwps[l, cell] + soil_w[l, cell] * soil_whcs[l, cell]) - soil_wsats[l, cell]
                soil_w[l, cell] -= grunoff / soil_whcs[l, cell]
                soil_lrunoff[l, cell] += grunoff
                lrunoff += grunoff
            end

            # percolation
            if (soil_w[l, cell] - percthres) > (1.0f-5 / soil_whcs[l, cell])
                # Calculate hydraulic conductivity
                HC = soil_Ks[l, cell] * ((soil_w[l, cell] * soil_whcs[l, cell] + soil_wpwps[l, cell]) / soil_wsats[l, cell])^soil_beta_soil[l, cell]
                # Calculate time constant
                TT = ((soil_w[l, cell] - percthres) * soil_whcs[l, cell]) / HC
                # Calculate percolation amount
                perc = ((soil_w[l, cell] - percthres) * soil_whcs[l, cell]) * (1 - exp(-24 / TT))
                # Correction of percolation for water content of the following layer
                if l < soil_layers
                    saturation_factor = 1 - (soil_w[l+1, cell] * soil_whcs[l+1, cell] + soil_w_fw[l+1, cell]) / (soil_wsats[l+1, cell] - soil_wpwps[l+1, cell])
                    if saturation_factor < 0
                        perc = zero(T)
                    else
                        perc *= sqrt(saturation_factor)
                    end
                else
                    saturation_factor = 1 - (soil_w[l, cell] * soil_whcs[l, cell] + soil_w_fw[l, cell]) / (soil_wsats[l, cell] - soil_wpwps[l, cell])
                    if saturation_factor < 0
                        perc = zero(T)
                    else
                        perc *= sqrt(saturation_factor)
                    end
                end
                
                soil_w[l, cell] -= perc / soil_whcs[l, cell] 

                if soil_w[l, cell] < 0
                    perc += soil_w[l, cell] * soil_whcs[l, cell]
                    soil_w[l, cell] = zero(T)
                end
                
                if l == soil_layers
                    soil_outflux_f[cell] += perc
                    soil_w_outflux[l, cell] += perc
                else
                    influx = perc
                    soil_w_influx[l+1, cell] += perc
                    soil_w_outflux[l, cell] += perc
                end

                soil_NO3[l, cell] += NO3perc_ly
                NO3perc_ly = zero(T)
                concNO3_mobile = zero(T)
                # determination of nitrate concentration in mobile water
                w_mobile = perc + srunoff + lrunoff
                if w_mobile > 1.0e-7
                    ww = -w_mobile / ((1 - anion_excl) * soil_wsats[l, cell])
                    vno3 = soil_NO3[l, cell] * (1 - exp(ww))
                    concNO3_mobile = max(vno3 / w_mobile, zero(0))
                end
                # NO3surf = zero(T)
                srunoff = zero(T)
                NO3lat = zero(T)
                if l == 1
                    NO3lat = NPERCO * concNO3_mobile * lrunoff
                else
                    NO3lat = concNO3_mobile * lrunoff
                end
                NO3lat = min(NO3lat, soil_NO3[l, cell])
                soil_NO3[l, cell] -= NO3lat

                # nitrate percolating from this layer
                NO3perc_ly = concNO3_mobile * perc
                NO3perc_ly = min(NO3perc_ly, soil_NO3[l, cell])
                soil_NO3[l, cell] -= NO3perc_ly
            end
        end
    end
    
    for l in 1:soil_layers
        soil_perc[l, cell] = soil_w_influx[l, cell] - soil_w_outflux[l, cell]
    end
    
end