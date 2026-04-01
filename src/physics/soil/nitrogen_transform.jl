function nitrogen_transform!(soil::Soil,
                             parm::LPJmLParam,
                             c_shift_fast::AbstractArray{T},
                             c_shift_slow::AbstractArray{T},
                             k_l = 0.0f0 # Parton et al., 2001 equ. 2
) where {T <: AbstractFloat}

    @unpack fastfrac, atmfrac, k_soil10 = parm

    # NO3 and N2O from mineralization of litter organic matter
    F_Nmineral = sum(soil.decom_litn, dims = 1) * atmfrac .* (fastfrac * c_shift_fast + (1.0f0 - fastfrac) * c_shift_slow);
    soil.NH4 .+= F_Nmineral * (1 - k_l)
    soil.NO3 .+= F_Nmineral * k_l

    # NO3 and N2O from mineralization of soil organic matter
    F_Nmineral = soil.decom_fastn + fsoil.decom_slown
    soil.NH4 .+= F_Nmineral * (1 - k_l)
    soil.NO3 .+= F_Nmineral * k_l

    # immobilization of N
    backend = KernelAbstractions.get_backend(soil.NH4)
    kernel = immobilize_kernel!(backend)
    kernel(soil.NH4, soil.NO3, c_shift_fast, c_shift_slow, soil.layer_depth, parm, ndrange=size(soil.NH4, 1))
    KernelAbstractions.synchronize(backend)

    # NO3 and N2O from nitrification
    backend = KernelAbstractions.get_backend(soil.NH4)
    kernel = nitrify_kernel!(backend)
    kernel(soil.NH4, soil.NO3, soil.swc, soil.wsats, soil.temp, soil.ph, parm, ndrange=size(soil.NH4, 1))
    KernelAbstractions.synchronize(backend)

end



@kernel function immobilize_kernel!(soil_NH4::AbstractArray{M},           
                                    soil_NO3::AbstractArray{M},
                                    c_shift_fast::AbstractArray{T},
                                    c_shift_slow::AbstractArray{T},
                                    soil_layer_depth::AbstractArray{T}, 
                                    parm::LPJmLParam;
                                    cn_ratio = 15,
                                    soil_layers = 5,
                                    k_N = 5f-3
) where {T <: AbstractFloat, M <: AbstractFloat}
    
    cell = @index(Global)

    @unpack fastfrac, atmfrac, k_soil10 = parm

    for l in 1:soil_layers

        N_sum = soil_NH4[l, cell] + soil_NO3[l, cell]
        if(N_sum > 0) # immobilization of N 
            n_immo = fastfrac * (1 - atmfrac) * (decom_sum_litc / cn_ratio - decom_sum_litn) * c_shift_fast[l] * N_sum / soil_layer_depth[l] * 1f3 / (k_N + N_sum / soil_layer_depth[l] * 1f3)
            if(n_immo > 0)
                if(n_immo > N_sum)
                    n_immo = N_sum
                end
                soil_slown[l, cell] += n_immo
                soil_NH4[l, cell] -= n_immo * soil_NH4[l, cell] / N_sum
                soil_NO3[l, cell] -= n_immo * soil_NO3[l, cell] / N_sum
            end
        end

        N_sum = soil_NH4[l, cell] + soil_NO3[l, cell]
        if(N_sum > 0) # immobilization of N 
            n_immo = (1 - fastfrac) * (1 - atmfrac) * (decom_sum_litc / cn_ratio - decom_sum_litn) * c_shift_slow[l] * N_sum / soil_layer_depth[l] * 1f3 / (k_N + N_sum / soil_layer_depth[l] * 1f3)
            if(n_immo > 0)
                if(n_immo > N_sum)
                    n_immo = N_sum
                end
                soil_slown[l, cell] += n_immo
                soil_NH4[l, cell] -= n_immo * soil_NH4[l, cell] / N_sum
                soil_NO3[l, cell] -= n_immo * soil_NO3[l, cell] / N_sum
            end
        end
    end

end


@kernel function nitrify_kernel!(soil_NH4::AbstractArray{M},           
                                 soil_NO3::AbstractArray{M},
                                 soil_swc::AbstractArray{M},
                                 soil_wsats::AbstractArray{M},
                                 soil_temp::AbstractArray{M},
                                 soil_ph::AbstractArray{T},
                                 parm::LPJmLParam;
                                 a_nit = 0.45f0,
                                 b_nit = 1.27f0,
                                 c_nit = 0.0012f0,
                                 d_nit = 2.84f0
) where {T <: AbstractFloat, M <: AbstractFloat}
    
    cell = @index(Global)

    @unpack k_max, k_2 = parm

    for l in 1:soil_layers

        x = soil_swc[l, cell] / soil_wsats[l, cell]
        n_nit = a_nit - b_nit
        m_nit = a_nit - c_nit
        z_nit = d_nit * (b_nit - a_nit) / (a_nit - c_nit)
        fac_wfps = ((x - b_nit) / n_nit)^(z_nit) * ((x - c_nit) / m_nit)^(d_nit)
        fac_temp = exp(-(soil_temp[l, cell] - T(18.79))^2 / T(2*5.26*5.26))

        F_NO3 = k_max * soil_NH4[l, cell] * fac_temp * fac_wfps * soil_ph[cell]
        if F_NO3 > soil_NH4[l, cell]
            F_NO3 = soil_NH4[l, cell]
        end
        # F_N2O = k_2 * F_NO3
        soil_NO3[l, cell] += F_NO3 * (1 - k_2)
        soil_NH4[l, cell] -= F_NO3
    end
end