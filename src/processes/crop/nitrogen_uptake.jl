"""
nuptake_crop!(crop, CFT, soil)

Compute root uptake of mineral nitrogen from soil NH4/NO3 pools. When
`auto_fertilizer=true`, supply the remaining plant demand as an explicit
external N input after root uptake, following LPJmL's AUTO_FERTILIZER mode.
"""
function nuptake_crop!(crop,
                       CFT::CFTParameters,
                       soil;
                       auto_fertilizer::Bool = false,
                       lpjmlparams::LPJmLParams = lpjmlparams
)

    kernel_params = (
        lpjmlparams = lpjmlparams,
        soil_layers = 5,
        auto_fertilizer = auto_fertilizer,
    )

    launch_1D!(
        nuptake_crop_kernel!,
        crop_prognostic(crop).nitrogen.total,
        crop_fluxes(crop).nitrogen.uptake,
        crop_fluxes(crop).nitrogen.auto_fertilizer,
        crop_prognostic(crop).nitrogen.leaf,
        crop_prognostic(crop).carbon.leaf,
        crop_prognostic(crop).nitrogen.root,
        crop_prognostic(crop).carbon.root,
        crop_stress_auxiliary(crop).nitrogen_demand_leaf,
        crop_stress_auxiliary(crop).nitrogen_demand_total,
        crop_prognostic(crop).nitrogen.sufficiency,
        crop_root_input(crop).distribution,
        crop_prognostic(crop).phenology.is_growing,
        soil_water_auxiliary(soil).relative_content,
        soil_water_prognostic(soil).saturation_fraction,
        soil_nitrogen_prognostic(soil).nitrate,
        soil_nitrogen_prognostic(soil).ammonium,
        soil_properties(soil).layer_depth,
        soil_thermal_prognostic(soil).temperature,
        kernel_constant(
            crop_prognostic(crop).nitrogen.total,
            (; CFT, kernel_params),
        ),
    )

end

"""Compute LPJmL's bounded soil-temperature response for root nitrogen uptake."""
@inline function compute_nitrogen_uptake_temperature_response(
    temperature::T,
    lower_temperature::T,
    optimum_temperature::T,
    reference_temperature::T,
) where {T <: AbstractFloat}
    return max(
        (temperature - lower_temperature) *
        (T(2) * optimum_temperature - lower_temperature - temperature) /
        (reference_temperature - lower_temperature) /
        (T(2) * optimum_temperature - lower_temperature - reference_temperature),
        zero(T),
    )
end

"""Compute the root-density and temperature scaling for one soil layer."""
@inline function compute_root_nitrogen_uptake_factor(
    temperature_response::T,
    plant_nitrogen_factor::T,
    root_carbon::T,
    root_fraction::T,
) where {T <: AbstractFloat}
    return temperature_response * plant_nitrogen_factor * root_carbon * root_fraction / T(1000)
end

"""Compute a layer's capped Michaelis–Menten NO₃ or NH₄ uptake potential."""
@inline function compute_mineral_nitrogen_uptake_potential(
    available::T,
    relative_water::T,
    saturation_fraction::T,
    layer_depth::T,
    root_factor::T,
    uptake_parameters,
) where {T <: AbstractFloat}
    available > zero(T) || return zero(T)
    water_scaler = relative_water > T(1e-7) ? one(T) : zero(T)
    saturation = available * water_scaler /
                 (available * water_scaler + T(uptake_parameters.Km) *
                  saturation_fraction * layer_depth / T(1000))
    potential = T(uptake_parameters.vmax) *
                (T(uptake_parameters.kmin) + saturation) * root_factor
    return min(potential, available)
end

@kernel inbounds = true function nuptake_crop_kernel!(
                                      crop_nitrogen::AbstractArray{T},
                                      crop_nuptake::AbstractArray{T},
                                      crop_nautofertilizer::AbstractArray{T},
                                      crop_leafn::AbstractArray{T},
                                      crop_leafc::AbstractArray{T},
                                      crop_rootn::AbstractArray{T},
                                      crop_rootc::AbstractArray{T},
                                      crop_ndemand_leaf::AbstractArray{T},
                                      crop_ndemand_tot::AbstractArray{T},
                                      crop_vscal::AbstractArray{T},
                                      crop_rootdist::AbstractArray{T},
                                      crop_isgrowing::AbstractArray{S},
                                      soil_w::AbstractArray{M},
                                      soil_wsat::AbstractArray{M},
                                      soil_NO3::AbstractArray{M},
                                      soil_NH4::AbstractArray{M},
                                      soil_layer_depth::AbstractArray{T},
                                      soil_temp::AbstractArray{M},
                                      fixed_parameters,
) where {T <: AbstractFloat, M <: AbstractFloat, S <: Integer}

    cell = @index(Global)
    parameters = kernel_value(fixed_parameters)
    CFT = parameters.CFT
    kernel_params = parameters.kernel_params

    @unpack lpjmlparams, soil_layers, auto_fertilizer = kernel_params

    @unpack T_0, T_m, T_r = lpjmlparams
    @unpack ncleaf, knstore, no3_uptake, nh4_uptake = CFT

    if crop_isgrowing[cell] == 1
        crop_nuptake[cell] = zero(T)
        crop_nautofertilizer[cell] = zero(T)

        mobile_carbon = crop_leafc[cell] + crop_rootc[cell]
        NCplant = mobile_carbon > zero(T) ?
                  (crop_leafn[cell] + crop_rootn[cell]) / mobile_carbon : T(ncleaf.low)
        nc_reference = T(2) / (one(T) / T(ncleaf.low) + one(T) / T(ncleaf.high))
        f_NCplant = clamp(
            (NCplant - T(ncleaf.high)) / (nc_reference - T(ncleaf.high)),
            zero(T),
            one(T),
        )

        leaf_nc = crop_leafc[cell] > zero(T) ?
                  crop_leafn[cell] / crop_leafc[cell] : zero(T)
        total_potential_uptake = zero(T)

        if leaf_nc < T(ncleaf.high) * (one(T) + T(knstore))
            # First pass: independent potential NO3 and NH4 uptake per layer.
            for l in 1:soil_layers
                temp_response = compute_nitrogen_uptake_temperature_response(
                    soil_temp[l, cell], T(T_0), T(T_m), T(T_r),
                )
                root_factor = compute_root_nitrogen_uptake_factor(
                    temp_response, f_NCplant, crop_rootc[cell], crop_rootdist[l],
                )
                total_potential_uptake += compute_mineral_nitrogen_uptake_potential(
                    max(zero(T), soil_NO3[l, cell]), soil_w[l, cell], soil_wsat[l, cell],
                    soil_layer_depth[l], root_factor, no3_uptake,
                )
                total_potential_uptake += compute_mineral_nitrogen_uptake_potential(
                    max(zero(T), soil_NH4[l, cell]), soil_w[l, cell], soil_wsat[l, cell],
                    soil_layer_depth[l], root_factor, nh4_uptake,
                )
            end
        end

        remaining_demand = max(zero(T), crop_ndemand_tot[cell] - crop_nitrogen[cell])
        n_uptake = min(total_potential_uptake, remaining_demand)

        if n_uptake > zero(T) && total_potential_uptake > zero(T)
            uptake_scale = n_uptake / total_potential_uptake

            # Second pass: remove exactly the accepted uptake from each pool.
            for l in 1:soil_layers
                temp_response = compute_nitrogen_uptake_temperature_response(
                    soil_temp[l, cell], T(T_0), T(T_m), T(T_r),
                )
                root_factor = compute_root_nitrogen_uptake_factor(
                    temp_response, f_NCplant, crop_rootc[cell], crop_rootdist[l],
                )
                no3_available = max(zero(T), soil_NO3[l, cell])
                if no3_available > zero(T)
                    no3_potential = compute_mineral_nitrogen_uptake_potential(
                        no3_available, soil_w[l, cell], soil_wsat[l, cell],
                        soil_layer_depth[l], root_factor, no3_uptake,
                    )
                    soil_NO3[l, cell] = max(
                        zero(T), soil_NO3[l, cell] - no3_potential * uptake_scale,
                    )
                end
                nh4_available = max(zero(T), soil_NH4[l, cell])
                if nh4_available > zero(T)
                    nh4_potential = compute_mineral_nitrogen_uptake_potential(
                        nh4_available, soil_w[l, cell], soil_wsat[l, cell],
                        soil_layer_depth[l], root_factor, nh4_uptake,
                    )
                    soil_NH4[l, cell] = max(
                        zero(T), soil_NH4[l, cell] - nh4_potential * uptake_scale,
                    )
                end
            end

            crop_nitrogen[cell] += n_uptake
            crop_nuptake[cell] = n_uptake
        end

        ndemand_leaf_opt = crop_ndemand_leaf[cell]
        nitrogen_deficit = max(zero(T), crop_ndemand_tot[cell] - crop_nitrogen[cell])
        if auto_fertilizer && nitrogen_deficit > zero(T)
            crop_nitrogen[cell] += nitrogen_deficit
            crop_nuptake[cell] += nitrogen_deficit
            crop_nautofertilizer[cell] = nitrogen_deficit
            crop_vscal[cell] = one(T)
        elseif nitrogen_deficit > zero(T)
            crop_ndemand_leaf[cell] = crop_leafn[cell]
            if ndemand_leaf_opt < T(1e-7)
                crop_vscal[cell] = one(T)
            else
                crop_vscal[cell] = min(
                    one(T),
                    crop_ndemand_leaf[cell] /
                    (ndemand_leaf_opt / (one(T) + T(knstore))),
                )
            end
        else
            crop_vscal[cell] = one(T)
        end

    else
        crop_nitrogen[cell] = zero(T)
        crop_nuptake[cell] = zero(T)
        crop_nautofertilizer[cell] = zero(T)
        # A deleted LPJmL stand has no active limitation. Keep the reusable
        # placeholder at the neutral multiplicative value.
        crop_vscal[cell] = one(T)
    end
end
