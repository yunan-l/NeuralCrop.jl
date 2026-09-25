"""
transpiration!(photos_adtmm, CFT, crop, pet, soil, co2; lpjmlparams=lpjmlparams)

Compute water demand/supply balance and layer-resolved transpiration uptake.
"""
function transpiration!(photos_adtmm::AbstractArray{T},
                        CFT::CFTParameters,
                        crop,
                        pet::PetPar,
                        soil,
                        co2::AbstractArray{T};
                        lpjmlparams::LPJmLParams = lpjmlparams
) where {T <: AbstractFloat}

    # Root-zone weighted water availability is accumulated inside the cell
    # kernel, avoiding a separate broadcast and reduction array every day.
    # supply = emax * wr .* (1 .- exp.(-0.04f0 * crop_prognostic(crop).carbon.root))
    # demand = ifelse.(crop_canopy_auxiliary(crop).canopy_conductance .> 0, (1 .- crop_canopy_auxiliary(crop).canopy_wet) .* pet.eeq * ALPHAM ./ (1 .+ (GM * ALPHAM) ./ crop_canopy_auxiliary(crop).canopy_conductance), zero(T))
    # transp = ifelse.(wr .> 0, min.(supply, demand) ./ wr .* fpc, zero(T)) # here the crop.fpc = 1, so we just omit it in the kernel fucntion

    kernel_params = (lpjmlparams = lpjmlparams, soil_layers = 5)

    launch_1D!(water_demand_supply_kernel!,
               crop_canopy_auxiliary(crop).canopy_conductance,
               photos_adtmm,
               co2,
               pet.daylength,
               crop_canopy_auxiliary(crop).fpar,
               crop_fluxes(crop).water.transpiration_layer,
               crop_prognostic(crop).water.demand_sum,
               crop_prognostic(crop).water.supply_sum,
               crop_stress_auxiliary(crop).water_deficit,
               crop_prognostic(crop).water.sufficiency,
               crop_prognostic(crop).carbon.root,
               crop_canopy_auxiliary(crop).canopy_wet,
               crop_prognostic(crop).phenology.is_growing,
               pet.eeq,
               crop_root_input(crop).distribution,
               crop_root_auxiliary(crop).zone_available_water,
               soil_water_auxiliary(soil).relative_content,
               soil_water_auxiliary(soil).holding_capacity_storage,
               kernel_constant(
                   crop_canopy_auxiliary(crop).canopy_conductance,
                   (; CFT, kernel_params),
               ))

end

"""Compute LPJmL canopy conductance from water-limited assimilation."""
@inline function compute_canopy_conductance(
    water_limited_assimilation::T,
    co2::T,
    daylength::T,
    fpar::T,
    minimum_conductance::T,
    lambda_optimum::T,
) where {T <: AbstractFloat}
    co2_bar = co2 * T(1e-5)
    co2_bar > zero(T) && daylength > zero(T) || return zero(T)
    denominator = co2_bar * (one(T) - lambda_optimum) * hour2sec(daylength)
    return T(1.6) * water_limited_assimilation / denominator +
           minimum_conductance * fpar
end

"""Compute root-water-limited transpiration supply for one crop column."""
@inline compute_transpiration_supply(
    maximum_supply::T, root_water::T, root_carbon::T,
) where {T <: AbstractFloat} =
    maximum_supply * root_water * (one(T) - exp(T(-0.04) * root_carbon))

"""Compute canopy transpiration demand from conductance and wet-canopy fraction."""
@inline function compute_transpiration_demand(
    canopy_wet::T,
    equilibrium_evaporation::T,
    alpha::T,
    conductance_shape::T,
    conductance::T,
) where {T <: AbstractFloat}
    conductance > zero(T) || return zero(T)
    return (one(T) - canopy_wet) * equilibrium_evaporation * alpha /
           (one(T) + (conductance_shape * alpha) / conductance)
end

"""Return the LPJmL 0–100 seasonal water-sufficiency diagnostic."""
@inline function compute_water_sufficiency(
    supplied::T, demanded::T,
) where {T <: AbstractFloat}
    demanded > zero(T) || return T(100)
    return clamp(T(100) * supplied / demanded, zero(T), T(100))
end

"""Cap one layer's transpiration extraction by its plant-available water."""
@inline function compute_layer_transpiration(
    transpiration::T,
    root_fraction::T,
    relative_water::T,
    holding_storage::T,
) where {T <: AbstractFloat}
    unconstrained = transpiration * root_fraction * relative_water
    capacity = relative_water * holding_storage
    return min(unconstrained, capacity), unconstrained > capacity
end

"""Recover actual canopy conductance after layer-wise uptake capping."""
@inline function compute_actual_canopy_conductance(
    current_conductance::T,
    actual_supply::T,
    demand::T,
    canopy_wet::T,
    equilibrium_evaporation::T,
    alpha::T,
    conductance_shape::T,
) where {T <: AbstractFloat}
    actual_supply < demand && equilibrium_evaporation > zero(T) || return current_conductance
    denominator = (one(T) - canopy_wet) * equilibrium_evaporation * alpha - actual_supply
    return denominator > zero(T) ? conductance_shape * alpha * actual_supply / denominator : zero(T)
end

@kernel inbounds = true function water_demand_supply_kernel!(
                                             crop_gp::AbstractArray{T},
                                             photos_adtmm::AbstractArray{T},
                                             co2::AbstractArray{T},
                                             daylength::AbstractArray{T},
                                             crop_fpar::AbstractArray{T},
                                             crop_trans_layer::AbstractArray{T},
                                             crop_w_demandsum::AbstractArray{T},
                                             crop_w_supplysum::AbstractArray{T},
                                             crop_wdf::AbstractArray{T},
                                             crop_wscal::AbstractArray{T},
                                             crop_rootc::AbstractArray{T},
                                             crop_canopy_wet::AbstractArray{T},
                                             crop_isgrowing::AbstractArray{S},
                                             pet_eeq::AbstractArray{T},
                                             crop_rootdist::AbstractArray{T},
                                             crop_rootzone_available_water::AbstractArray{T},
                                             soil_w::AbstractArray{M},
                                             soil_whcs::AbstractArray{M},
                                             fixed_parameters,
) where {T <: AbstractFloat, M <: AbstractFloat, S <: Integer}

    cell = @index(Global)
    parameters = kernel_value(fixed_parameters)
    CFT = parameters.CFT
    kernel_params = parameters.kernel_params

    @unpack lpjmlparams, soil_layers = kernel_params
    @unpack ALPHAM, GM, LAMBDA_OPT = lpjmlparams
    @unpack fpc, emax, gmin = CFT

    co2_index = length(co2) == 1 ? 1 : cell
    crop_gp[cell] = compute_canopy_conductance(
        photos_adtmm[cell], co2[co2_index], daylength[cell], crop_fpar[cell],
        T(gmin), T(LAMBDA_OPT),
    )

    wr = zero(T)
    rootzone_water = zero(T)
    for l in 1:soil_layers
        wr += soil_w[l, cell] * crop_rootdist[l]
        if l <= 3
            rootzone_water += soil_w[l, cell] * soil_whcs[l, cell] * crop_rootdist[l]
        end
    end
    crop_rootzone_available_water[cell] = rootzone_water

    if crop_isgrowing[cell] == 1
        supply = compute_transpiration_supply(T(emax), wr, crop_rootc[cell])
        demand = compute_transpiration_demand(
            crop_canopy_wet[cell], pet_eeq[cell], T(ALPHAM), T(GM), crop_gp[cell],
        )

        crop_w_demandsum[cell] += demand
        if supply > demand
            crop_w_supplysum[cell] += demand
        else
            crop_w_supplysum[cell] += supply
        end

        crop_wdf[cell] = compute_water_sufficiency(
            crop_w_supplysum[cell], crop_w_demandsum[cell],
        )

        if pet_eeq[cell] > 0.0 && crop_gp[cell] > 0.0
            crop_wscal[cell] = (emax * wr) / (pet_eeq[cell] * ALPHAM / (one(T) + (GM * ALPHAM) / crop_gp[cell]))
            if crop_wscal[cell] > 1.0
                crop_wscal[cell] = one(T)
            end
        else
            crop_wscal[cell] = one(T)
        end

        # Potential transpiration constrained by demand/supply and canopy fraction.
        if wr > 0
            transp = min(supply, demand) / wr * fpc
        else
            transp = zero(T)
        end

        transp_cor = zero(T)

        # Apply layer-wise extraction cap so uptake does not exceed layer storage.
        if transp > 0
            for l in 1:soil_layers
                transp_tmp, capped = compute_layer_transpiration(
                    transp, crop_rootdist[l], soil_w[l, cell], soil_whcs[l, cell],
                )
                transp_cor += transp_tmp
                capped && transp_cor < T(1e-5) && (transp_cor = zero(T))
            end
        else
            transp_cor = zero(T)
        end

        if wr > 0
            transp = transp_cor / wr
        else
            transp = zero(T)
        end

        # LPJmL recomputes actual canopy conductance after layer extraction.
        # Store it in `gp`; downstream lambda solving consumes this actual value.
        actual_supply = fpc > zero(T) ? transp_cor / fpc : zero(T)
        crop_gp[cell] = compute_actual_canopy_conductance(
            crop_gp[cell], actual_supply, demand, crop_canopy_wet[cell], pet_eeq[cell],
            T(ALPHAM), T(GM),
        )

        # Distribute corrected transpiration back to layers by root distribution.
        for l in 1:soil_layers
            crop_trans_layer[l, cell], _ = compute_layer_transpiration(
                transp, crop_rootdist[l], soil_w[l, cell], soil_whcs[l, cell],
            )
        end
    else
        crop_gp[cell] = zero(T)
        for l in 1:soil_layers
            crop_trans_layer[l, cell] = zero(T)
        end
        crop_w_demandsum[cell] = zero(T)
        crop_w_supplysum[cell] = zero(T)
        crop_wdf[cell] = zero(T)
        # Neutral stress for an absent stand; is_growing still gates all fluxes.
        crop_wscal[cell] = one(T)
    end
end
