"""
evaporation!(pet_eeq, crop, soil)

Compute layer-wise bare-soil evaporation constrained by near-surface water.
"""
function evaporation!(pet_eeq::AbstractArray{T},
                      crop,
                      soil;
                      lpjmlparams::LPJmLParams = lpjmlparams

) where {T <: AbstractFloat}

    kernel_params = (lpjmlparams = lpjmlparams, soil_layers = 5)

    launch_1D!(evaporation_kernel!,
               pet_eeq,
               crop_canopy_auxiliary(crop).fpar,
               crop_fluxes(crop).water.transpiration_layer,
               crop_canopy_auxiliary(crop).canopy_wet,
               soil_water_auxiliary(soil).relative_content,
               soil_water_auxiliary(soil).free_water,
               soil_water_auxiliary(soil).holding_capacity_storage,
               soil_water_fluxes(soil).evaporation,
               soil_surface_litter_prognostic(soil).cover,
               soil_surface_litter_auxiliary(soil).water_capacity,
               soil_surface_litter_prognostic(soil).water_storage,
               soil_surface_litter_fluxes(soil).evaporation,
               soil_properties(soil).layer_depth,
               kernel_constant(pet_eeq, kernel_params))

end

"""
    compute_litter_evaporation(evaporation_energy, available_energy, storage,
                               capacity, cover)

Return the bounded evaporation loss from the surface-litter water store. The
wetness-squared dependence and all three physical caps follow the existing
LPJmL-compatible implementation.
"""
@inline function compute_litter_evaporation(evaporation_energy::A,
                                             available_energy::B,
                                             storage::C,
                                             capacity::D,
                                             cover::E) where {A <: AbstractFloat,
                                                               B <: AbstractFloat,
                                                               C <: AbstractFloat,
                                                               D <: AbstractFloat,
                                                               E <: AbstractFloat}
    # `evaporation_energy` can be Float64 because LPJmL's scalar
    # Priestley--Taylor parameter is retained at that precision. Do not narrow
    # it here: the original kernel intentionally used Julia promotion.
    (capacity > eps(D) && available_energy > zero(B)) || return zero(C)
    wetness = clamp(storage / capacity, zero(C), one(C))
    return min(evaporation_energy * wetness * wetness * cover, storage, available_energy)
end

"""
    compute_soil_evaporation_ratio(evaporation_energy, available_energy,
                                   available_water, holding_storage, litter_cover)

Compute the scalar fraction of near-surface liquid water removed by bare-soil
evaporation. Layer aggregation and layer-wise application remain in the kernel
to preserve the existing five-layer mass update.
"""
@inline function compute_soil_evaporation_ratio(evaporation_energy::A,
                                                 available_energy::B,
                                                 available_water::C,
                                                 holding_storage::D,
                                                 litter_cover::E) where {A <: AbstractFloat,
                                                                           B <: AbstractFloat,
                                                                           C <: AbstractFloat,
                                                                           D <: AbstractFloat,
                                                                           E <: AbstractFloat}
    (evaporation_energy > 1.0f-5 && available_energy > 1.0f-5 &&
     available_water > zero(C)) || return zero(C)
    soil_evaporation = evaporation_energy /
        (1 + exp(5 - 10 * available_water / holding_storage)) *
        max(0.05, 1 - litter_cover)
    return soil_evaporation / available_water
end

@kernel inbounds = true function evaporation_kernel!(
                                     pet_eeq::AbstractArray{T},
                                     crop_fpar::AbstractArray{T},
                                     crop_trans_layer::AbstractArray{M},
                                     crop_canopy_wet::AbstractArray{T},
                                     soil_w::AbstractArray{M},
                                     soil_w_fw::AbstractArray{M},
                                     soil_whcs::AbstractArray{M},
                                     soil_evap::AbstractArray{M},
                                     soil_agtop_cover::AbstractArray{T},
                                     litter_water_capacity::AbstractArray{T},
                                     litter_water_storage::AbstractArray{T},
                                     litter_evaporation::AbstractArray{T},
                                     soil_layer_depth::AbstractArray{T},
                                     fixed_parameters,
) where {T <: AbstractFloat, M <: AbstractFloat}

    cell = @index(Global)
    kernel_params = kernel_value(fixed_parameters)

    @unpack lpjmlparams, soil_layers = kernel_params

    @unpack PRIESTLEY_TAYLOR = lpjmlparams  # Priestley-Taylor coefficient

    soildepth_evap = lpjmlparams.soildepth_evap

    evap_energy = pet_eeq[cell] * PRIESTLEY_TAYLOR * max(1 - crop_fpar[cell], 0.05)
    # evap_litter = pet_eeq[cell] * PRIESTLEY_TAYLOR * (1 - crop_canopy_wet[cell]) - sum(crop_trans_layer[:, cell])

    crop_trans_layer_sum = zero(T)
    for l in 1:soil_layers
        crop_trans_layer_sum += crop_trans_layer[l, cell]
    end

    available_evaporation = max(
        pet_eeq[cell] * PRIESTLEY_TAYLOR *
        (one(T) - crop_canopy_wet[cell]) - crop_trans_layer_sum,
        zero(T),
    )
    litter_evaporation[cell] = compute_litter_evaporation(
        evap_energy, available_evaporation, litter_water_storage[cell],
        litter_water_capacity[cell], soil_agtop_cover[cell],
    )
    litter_water_storage[cell] -= litter_evaporation[cell]

    evap_ratio = zero(T)
    if evap_energy > 1.0f-5 && (pet_eeq[cell] * PRIESTLEY_TAYLOR * (1 - crop_canopy_wet[cell]) - crop_trans_layer_sum) > 1.0f-5
        # w_evap is water content in soildepth_evap that can evaporate
        w_evap = zero(T)
        whcs_evap = zero(T)

        for l in 1:soil_layers
            if soildepth_evap > 0
                fraction = min(1, soildepth_evap / soil_layer_depth[l])
                liquid_above_pwp = soil_w[l, cell] * soil_whcs[l, cell] +
                                   soil_w_fw[l, cell]
                w_evap += max(liquid_above_pwp - crop_trans_layer[l, cell], zero(T)) * fraction
                whcs_evap += soil_whcs[l, cell] * fraction
                soildepth_evap -= soil_layer_depth[l]
            end
        end

        evap_ratio = compute_soil_evaporation_ratio(
            evap_energy,
            pet_eeq[cell] * PRIESTLEY_TAYLOR * (one(T) - crop_canopy_wet[cell]) -
                crop_trans_layer_sum,
            w_evap, whcs_evap, soil_agtop_cover[cell],
        )
    end

    soildepth_evap = lpjmlparams.soildepth_evap
    for l in 1:soil_layers
        if soildepth_evap > 0
            fraction = min(1, soildepth_evap / soil_layer_depth[l])
            liquid_above_pwp = soil_w[l, cell] * soil_whcs[l, cell] +
                               soil_w_fw[l, cell]
            soil_evap[l, cell] = liquid_above_pwp * evap_ratio * fraction
            soildepth_evap -= soil_layer_depth[l]
        end
    end

end
