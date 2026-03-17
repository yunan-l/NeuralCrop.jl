function phenology_crop!(crop::Crop,
                         climbuf_V_req::AbstractArray{T},
                         PFT::PftParameters,
                         temp::AbstractArray{T},
                         daylength::AbstractArray{T},
) where {T <: AbstractFloat}
    
    @unpack laimax = PFT

    crop.harvesting0 .= crop.harvesting

    backend = get_backend(temp)
    
    kernel = phenology_kernel!(backend)
    
    kernel(crop.phu, 
           crop.vdsum, 
           crop.husum, 
           crop.fphu,
           crop.flaimax,
           crop.senescence, 
           crop.senescence0, 
           crop.harvesting, 
           crop.growingdays,
           crop.isgrowing,
           crop.wtype,
           PFT, 
           climbuf_V_req, 
           temp, 
           daylength, 
           ndrange=length(climbuf_V_req)
    )
    
    synchronize(backend)

    lai_crop!(crop, PFT)
    
    crop.phen = crop.lai / laimax
  
end


@kernel function phenology_kernel!(crop_phu::AbstractArray{T},           
                                   crop_vdsum::AbstractArray{T},           
                                   crop_husum::AbstractArray{T},            
                                   crop_fphu::AbstractArray{T}, 
                                   crop_flaimax::AbstractArray{T}, 
                                   crop_senescence::AbstractArray{B},
                                   crop_senescence0::AbstractArray{B},                  
                                   crop_harvesting::AbstractArray{B},         
                                   crop_growingdays::AbstractArray{S},
                                   crop_isgrowing::AbstractArray{S},
                                   crop_wtype::AbstractArray{B},
                                   PFT::PftParameters,
                                   climbuf_V_req::AbstractArray{T},
                                   temp::AbstractArray{T},
                                   daylength::AbstractArray{T},
) where {T <: AbstractFloat, B <: Bool, S <: Integer}
    
    cell = @index(Global)

    @unpack basetemp, tv_eff, tv_opt, fphuc, flaimaxc, fphuk, flaimaxk, fphusen, flaimaxharvest, psens, pb, ps, hlimit, sla, shapesenescencenorm = PFT

    crop_senescence0[cell] = crop_senescence[cell]
    
    if crop_isgrowing[cell] == 1
        crop_growingdays[cell] += 1
        if crop_husum[cell] < crop_phu[cell]
            hu = max(0, temp[cell] - basetemp.low)
            if crop_wtype[cell] # winter crops with vernalization requirements
                if crop_vdsum[cell] < climbuf_V_req[cell]
                    if temp[cell] >= tv_eff.low && temp[cell] < tv_opt.low
                        vd_inc = (temp[cell] - tv_eff.low) / (tv_opt.low - tv_eff.low)
                    elseif temp[cell] <= tv_eff.high && temp[cell] >= tv_opt.high
                        vd_inc = (tv_eff.high - temp[cell]) / (tv_eff.high - tv_opt.high)
                    elseif temp[cell] >= tv_opt.low && temp[cell] < tv_opt.high
                        vd_inc = one(T)
                    else
                        vd_inc = zero(T)
                    end
                else
                    vd_inc = zero(T)
                end
                
                crop_vdsum[cell] += max(zero(T), vd_inc)
                #Calculation of vernalization reduction factor
                vd_b = climbuf_V_req[cell] / 5 # base requirements, 20% of total vernalization requirements
                
                if crop_vdsum[cell] < vd_b
                    vrf = zero(T)
                elseif crop_vdsum[cell] >= vd_b && crop_vdsum[cell] < climbuf_V_req[cell]
                    vrf = max(zero(T), min(one(T), (crop_vdsum[cell] - vd_b) / (climbuf_V_req[cell] - vd_b)))
                else
                    vrf = one(T)
                end
                
            else
                vrf = one(T)
            end
            
            #Response to photoperiodism (still inactive, yet. This means that PFT.psens == 1 for all crops)
            if crop_fphu[cell] <= fphusen
                prf = (one(T) - psens) * min(one(T), max(zero(T), (daylength[cell] - pb) / (ps - pb))) + psens
            else
                prf = one(T)
            end
            
            #Calculation of temperature sum (deg Cd)
            crop_husum[cell] += hu * vrf * prf

            #fraction of growing season
            crop_fphu[cell] = min(one(T), crop_husum[cell] / crop_phu[cell])

            if crop_fphu[cell] < fphusen
                c = fphuc / flaimaxc - fphuc
                k = fphuk / flaimaxk - fphuk
                crop_flaimax[cell] = crop_fphu[cell] / (crop_fphu[cell] + c * (c/k) ^ ((fphuc - crop_fphu[cell]) / (fphuk - fphuc)))
            else 
                crop_senescence[cell] = true
                crop_flaimax[cell] = ((1 - crop_fphu[cell]) / (1 - fphusen)) ^ shapesenescencenorm * (1 - flaimaxharvest) + flaimaxharvest
            end

        else
            crop_harvesting[cell] = true
        end
        
        if(crop_growingdays[cell] == hlimit)
            crop_harvesting[cell] = true
        end
    else
        crop_vdsum[cell] = zero(T)
        crop_husum[cell] = zero(T)
        crop_fphu[cell] = zero(T)
        crop_senescence[cell] = false
        crop_growingdays[cell] = 0
        crop_flaimax[cell] = zero(T)
    end  
end