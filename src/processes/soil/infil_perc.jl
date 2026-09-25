"""
infil_perc!(soil; lpjmlparams=lpjmlparams)

Update infiltration, percolation, runoff, and nitrate transport through soil layers.
"""
function infil_perc!(soil;
                     lpjmlparams::LPJmLParams = lpjmlparams)
    return infil_perc!(
        soil,
        soil_water_fluxes(soil).infiltration,
        soil_snow_fluxes(soil).melt,
        soil_water_fluxes(soil).infiltration;
        lpjmlparams = lpjmlparams,
        transfer_heat = false,
    )
end

"""
    apply_percolation_enthalpy!(soil)

Apply the layer-wise `perc_energy` ledger after the corresponding water-stock
update. Thermal properties and phase fractions are rebuilt from the new water
mass, matching LPJmL's `apply_perc_enthalpy`/thermal-property reconciliation.
"""
function apply_percolation_enthalpy!(
    soil;
    thermalparams::SoilThermalParams = soil_thermal_params,
)
    launch_1D!(
        apply_percolation_enthalpy_kernel!,
        soil_thermal_fluxes(soil).percolation_energy_residual,
        soil_thermal_prognostic(soil).temperature,
        soil_thermal_prognostic(soil).enthalpy,
        soil_thermal_prognostic(soil).frozen_fraction,
        soil_thermal_prognostic(soil).freeze_depth,
        soil_thermal_prognostic(soil).heat_capacity_frozen,
        soil_thermal_prognostic(soil).heat_capacity_unfrozen,
        soil_thermal_prognostic(soil).latent_heat,
        soil_thermal_prognostic(soil).conductivity_frozen,
        soil_thermal_prognostic(soil).conductivity_unfrozen,
        soil_thermal_prognostic(soil).water_reference,
        soil_thermal_fluxes(soil).percolation_energy,
        soil_thermal_fluxes(soil).rain_energy_input,
        soil_thermal_fluxes(soil).snowmelt_energy_input,
        soil_thermal_fluxes(soil).lateral_runoff_energy_output,
        soil_thermal_fluxes(soil).bottom_drainage_energy_output,
        soil_thermal_input(soil).diffusivity_0,
        soil_water_prognostic(soil).storage,
        soil_water_prognostic(soil).ice_storage,
        soil_water_auxiliary(soil).saturation_storage,
        soil_properties(soil).layer_depth,
        kernel_constant(
            soil_thermal_fluxes(soil).percolation_energy_residual, thermalparams,
        ),
    )
    partition_soil_water_ice!(soil)
    return nothing
end

@kernel inbounds = true function apply_percolation_enthalpy_kernel!(
    percolation_energy_residual::AbstractArray{T},
    temperature::AbstractArray{T},
    enthalpy::AbstractArray{T},
    frozen_fraction::AbstractArray{T},
    freeze_depth::AbstractArray{T},
    heat_capacity_frozen::AbstractArray{T},
    heat_capacity_unfrozen::AbstractArray{T},
    latent_heat::AbstractArray{T},
    conductivity_frozen::AbstractArray{T},
    conductivity_unfrozen::AbstractArray{T},
    water_reference::AbstractArray{T},
    percolation_energy::AbstractArray{T},
    rain_energy_input::AbstractArray{T},
    snowmelt_energy_input::AbstractArray{T},
    lateral_energy_output::AbstractArray{T},
    bottom_energy_output::AbstractArray{T},
    diffusivity_0::AbstractArray{T},
    liquid_water::AbstractArray{T},
    ice_water::AbstractArray{T},
    saturation_storage::AbstractArray{T},
    layer_depth::AbstractArray{T},
    fixed_parameters,
) where {T <: AbstractFloat}
    cell = @index(Global)
    thermalparams = kernel_value(fixed_parameters)
    transfer_sum = zero(T)

    for layer in 1:5
        depth_m = max(layer_depth[layer] * T(0.001), eps(T))
        layer_transfer = percolation_energy[layer, cell]
        transfer_sum += layer_transfer
        updated_enthalpy = enthalpy[layer, cell] + layer_transfer / depth_m
        total_water = max(liquid_water[layer, cell] + ice_water[layer, cell], zero(T))
        cf, cu, lh, kf, ku = layer_thermal_properties(
            total_water,
            saturation_storage[layer, cell],
            layer_depth[layer],
            diffusivity_0[cell],
            thermalparams,
        )
        layer_temperature = enthalpy_temperature(updated_enthalpy, cf, cu, lh)
        layer_frozen_fraction = enthalpy_frozen_fraction(updated_enthalpy, lh)

        enthalpy[layer, cell] = updated_enthalpy
        temperature[layer, cell] = layer_temperature
        frozen_fraction[layer, cell] = layer_frozen_fraction
        freeze_depth[layer, cell] = layer_frozen_fraction * layer_depth[layer]
        heat_capacity_frozen[layer, cell] = cf
        heat_capacity_unfrozen[layer, cell] = cu
        latent_heat[layer, cell] = lh
        conductivity_frozen[layer, cell] = kf
        conductivity_unfrozen[layer, cell] = ku
        water_reference[layer, cell] = total_water
        ice_water[layer, cell] = layer_frozen_fraction * total_water
        liquid_water[layer, cell] = total_water - ice_water[layer, cell]
        percolation_energy[layer, cell] = zero(T)
    end

    boundary_energy = rain_energy_input[cell] + snowmelt_energy_input[cell] -
                      lateral_energy_output[cell] - bottom_energy_output[cell]
    percolation_energy_residual[cell] = transfer_sum - boundary_energy
end

"""
    infil_perc!(soil, precipitation, snowmelt, air_temperature;
                transfer_heat=true)

Route water, nitrate, and the enthalpy carried by liquid water. The upper
boundary follows LPJmL: rainfall enters at air temperature and meltwater at
0 °C. Water leaving a layer carries that layer's liquid-water enthalpy.
"""
function infil_perc!(soil,
                     precipitation::AbstractArray{T},
                     snowmelt::AbstractArray{T},
                     air_temperature::AbstractArray{T};
                     lpjmlparams::LPJmLParams = lpjmlparams,
                     thermalparams::SoilThermalParams{T} = soil_thermal_params,
                     transfer_heat::Bool = true) where {T <: AbstractFloat}
    # One-cell kernel launch; each thread updates the full vertical soil column for that cell.
    kernel_params = (;
        lpjmlparams,
        thermalparams,
        soil_layers = 5,
        anion_excl = T(0.3),
        NPERCO = T(0.4),
        transfer_heat,
    )

    water_fluxes = soil_water_fluxes(soil)
    water = soil_water_prognostic(soil)
    water_auxiliary = soil_water_auxiliary(soil)
    thermal = soil_thermal_prognostic(soil)
    thermal_fluxes = soil_thermal_fluxes(soil)
    launch_1D!(
        infil_perc_kernel!,
        water_fluxes.infiltration,
        (
            relative_content = water_auxiliary.relative_content,
            holding_capacity_storage = water_auxiliary.holding_capacity_storage,
            free_water = water_auxiliary.free_water,
            ice_storage = water.ice_storage,
            available_ice_storage = water.available_ice_storage,
            free_ice_storage = water.free_ice_storage,
            saturation_fraction = water.saturation_fraction,
            saturation_storage = water_auxiliary.saturation_storage,
            wilting_fraction = water_auxiliary.wilting_fraction,
            wilting_storage = water_auxiliary.wilting_storage,
        ),
        (
            influx = water_fluxes.influx,
            outflux = water_fluxes.outflux,
            saturated_conductivity = water_auxiliary.saturated_conductivity,
            surface_runoff = water_fluxes.surface_runoff,
            lateral_runoff = water_fluxes.lateral_runoff,
            bottom_drainage = water_fluxes.bottom_drainage,
            percolation = water_fluxes.percolation,
        ),
        (
            freeze_depth = thermal.freeze_depth,
            surface_cover = soil_surface_litter_prognostic(soil).cover,
            beta = water_auxiliary.beta,
            nitrate = soil_nitrogen_prognostic(soil).nitrate,
            nitrogen_leaching = soil_nitrogen_fluxes(soil).leaching,
            layer_depth = soil_properties(soil).layer_depth,
            sand_fraction = soil_properties(soil).sand_fraction,
            tillage_density_factor = soil_management_prognostic(soil).tillage_density_factor,
            temperature = thermal.temperature,
        ),
        (
            percolation = thermal_fluxes.percolation_energy,
            rain_input = thermal_fluxes.rain_energy_input,
            snowmelt_input = thermal_fluxes.snowmelt_energy_input,
            lateral_output = thermal_fluxes.lateral_runoff_energy_output,
            bottom_output = thermal_fluxes.bottom_drainage_energy_output,
        ),
        (; precipitation, snowmelt, air_temperature),
        kernel_constant(water_fluxes.infiltration, kernel_params),
    )
    return nothing
end

"""
    compute_infiltration_slug(slug, relative_storage, available_storage, exponent)

LPJmL's bounded surface-infiltration response for one small precipitation
slug. A saturated surface has no hydraulic intake; the remainder becomes
surface runoff in the calling kernel.
"""
@inline function compute_infiltration_slug(slug::T,
                                            relative_storage::T,
                                            available_storage::T,
                                            exponent::T) where {T <: AbstractFloat}
    saturation_deficit = one(T) - relative_storage / available_storage
    saturation_deficit >= zero(T) || return zero(T)
    return slug * saturation_deficit^(one(T) / exponent)
end

"""
    compute_liquid_water_enthalpy(temperature, water_heat_capacity,
                                  ice_heat_capacity, fusion_heat)

Return volumetric enthalpy transported by water leaving a soil layer. Liquid
water above freezing carries fusion heat plus sensible heat; frozen material
uses the ice heat capacity, matching the existing percolation-energy ledger.
"""
@inline function compute_liquid_water_enthalpy(temperature::T,
                                               water_heat_capacity::T,
                                               ice_heat_capacity::T,
                                               fusion_heat::T) where {T <: AbstractFloat}
    return temperature >= zero(T) ?
        fusion_heat + water_heat_capacity * temperature :
        ice_heat_capacity * temperature
end

"""
    compute_percolation(storage, holding_capacity, threshold, conductivity,
                         saturation_fraction, beta, ice_fraction,
                         receiver_saturation_factor)

Compute a daily LPJmL hydraulic drainage amount before the donor-stock clamp.
The receiver factor is supplied by the caller because it depends on the next
layer (or the bottom layer for free drainage).
"""
@inline function compute_percolation(storage::T,
                                     holding_capacity::T,
                                     threshold::T,
                                     conductivity::T,
                                     saturation_fraction::T,
                                     beta::T,
                                     ice_fraction::T,
                                     receiver_saturation_factor::T) where {T <: AbstractFloat}
    excess = (storage - threshold) * holding_capacity
    excess > T(1e-5) || return zero(T)
    hydraulic_conductivity = conductivity * saturation_fraction^beta *
        T(10)^(-T(6) * ice_fraction)
    time_constant = excess / hydraulic_conductivity
    potential = excess * (one(T) - exp(-T(24) / time_constant))
    return receiver_saturation_factor < zero(T) ? zero(T) :
        potential * sqrt(receiver_saturation_factor)
end

"""
    compute_mobile_nitrate_concentration(nitrate, mobile_water, exclusion,
                                         saturation_storage)

Calculate nitrate concentration in the mobile-water fraction. The threshold
and non-negative clamp preserve the previous LPJmL-compatible behavior when a
dry layer has no mobile flux.
"""
@inline function compute_mobile_nitrate_concentration(nitrate::T,
                                                       mobile_water::T,
                                                       anion_exclusion::T,
                                                       saturation_storage::T) where {T <: AbstractFloat}
    mobile_water > T(1e-7) || return zero(T)
    exponent = -mobile_water / ((one(T) - anion_exclusion) * saturation_storage)
    dissolved = nitrate * (one(T) - exp(exponent))
    return max(dissolved / mobile_water, zero(T))
end

"""
    compute_tillage_settling_fraction(infiltration, sand_fraction, depth)

Return the rainfall-driven approach fraction toward untilled bulk density for
the top layer. This is a state-independent LPJmL scalar response; its actual
state update remains in the kernel after today's infiltration has completed.
"""
@inline function compute_tillage_settling_fraction(infiltration::T,
                                                    sand_fraction::T,
                                                    depth::T) where {T <: AbstractFloat}
    sand_percent = sand_fraction * T(100)
    settling_index = T(0.2) * infiltration * (
        one(T) + T(2) * sand_percent /
        (sand_percent + exp(T(8.597) - T(0.075) * sand_percent))
    ) / max(depth * T(0.001), eps(T))^T(0.6)
    return settling_index / (settling_index + exp(T(3.92) - T(0.0226) * settling_index))
end

@kernel inbounds = true function infil_perc_kernel!(
                                    infil::AbstractArray{T},
                                    water,
                                    water_fluxes,
                                    soil,
                                    energy,
                                    forcing,
                                    fixed_parameters,
) where {T <: AbstractFloat}

    cell = @index(Global)
    kernel_params = kernel_value(fixed_parameters)

    soil_w = water.relative_content
    soil_whcs = water.holding_capacity_storage
    soil_w_fw = water.free_water
    soil_ice = water.ice_storage
    soil_ice_available = water.available_ice_storage
    soil_ice_free = water.free_ice_storage
    soil_wsat = water.saturation_fraction
    soil_wsats = water.saturation_storage
    soil_wpwp = water.wilting_fraction
    soil_wpwps = water.wilting_storage
    soil_w_influx = water_fluxes.influx
    soil_w_outflux = water_fluxes.outflux
    soil_Ks = water_fluxes.saturated_conductivity
    soil_srunoff = water_fluxes.surface_runoff
    soil_lrunoff = water_fluxes.lateral_runoff
    soil_outflux_f = water_fluxes.bottom_drainage
    soil_perc = water_fluxes.percolation
    soil_freeze_depth = soil.freeze_depth
    soil_agtop_cover = soil.surface_cover
    soil_beta_soil = soil.beta
    soil_NO3 = soil.nitrate
    soil_n_leaching = soil.nitrogen_leaching
    soil_layer_depth = soil.layer_depth
    soil_sand_fraction = soil.sand_fraction
    tillage_density_factor = soil.tillage_density_factor
    soil_temperature = soil.temperature
    soil_perc_energy = energy.percolation
    rain_energy_input = energy.rain_input
    snowmelt_energy_input = energy.snowmelt_input
    lateral_energy_output = energy.lateral_output
    bottom_energy_output = energy.bottom_output
    precipitation = forcing.precipitation
    snowmelt = forcing.snowmelt
    air_temperature = forcing.air_temperature

    @unpack lpjmlparams, thermalparams, soil_layers, anion_excl, NPERCO,
            transfer_heat = kernel_params

    @unpack soil_infil, soil_infil_litter, percthres = lpjmlparams

    freewater = zero(T)
    soil_srunoff[cell] = zero(T)
    soil_outflux_f[cell] = zero(T)
    soil_n_leaching[cell] = zero(T)
    rain_energy_input[cell] = zero(T)
    snowmelt_energy_input[cell] = zero(T)
    lateral_energy_output[cell] = zero(T)
    bottom_energy_output[cell] = zero(T)
    for l in 1:soil_layers
        freewater += soil_w_fw[l, cell]
        if soil_w[l, cell] + soil_ice_available[l, cell] / soil_whcs[l, cell] > one(T)
            freewater += (soil_w[l, cell] +
                          soil_ice_available[l, cell] / soil_whcs[l, cell] - one(T)) *
                         soil_whcs[l, cell]
        end
        soil_lrunoff[l, cell] = zero(T)
        soil_w_influx[l, cell] = zero(T)
        soil_w_outflux[l, cell] = zero(T)
        soil_perc_energy[l, cell] = zero(T)
    end

    @unpack water_heat_capacity, ice_heat_capacity,
            volumetric_fusion_heat = thermalparams
    total_top_water = max(precipitation[cell], zero(T))
    melt_top_water = min(max(snowmelt[cell], zero(T)), total_top_water)
    rain_top_water = max(total_top_water - melt_top_water, zero(T))
    top_water_denominator = rain_top_water + melt_top_water
    # Keep inactive dry-day divisions finite for reverse AD. Wet days retain
    # the original denominator and threshold.
    safe_top_water_denominator = max(top_water_denominator, eps(T))
    rain_fraction = top_water_denominator > eps(T) ?
                    rain_top_water / safe_top_water_denominator : zero(T)
    melt_fraction = top_water_denominator > eps(T) ?
                    melt_top_water / safe_top_water_denominator : zero(T)
    rain_volumetric_enthalpy = volumetric_fusion_heat +
                               water_heat_capacity * air_temperature[cell]
    melt_volumetric_enthalpy = volumetric_fusion_heat
    top_volumetric_enthalpy = rain_fraction * rain_volumetric_enthalpy +
                              melt_fraction * melt_volumetric_enthalpy

    soil_infil *= (one(T) + soil_agtop_cover[cell] * soil_infil_litter)
    influx = zero(T)
    iter = 0
    # LPJmL bounds the slug loop at MAXITER = 1000. This is a safety cap,
    # not a convergence criterion: any unprocessed rainfall becomes surface
    # runoff when the cap is reached.
    while (infil[cell] > T(1e-5) || freewater > T(1e-5)) && iter < 1000
        iter += 1
        NO3perc_ly = zero(T)
        freewater = zero(T)
        # Process infiltration in bounded slugs for numerical stability.
        slug = min(4, infil[cell])
        infil[cell] -= slug

        # Calculate influx to first soil layer
        top_storage = soil_w[1, cell] * soil_whcs[1, cell] + soil_w_fw[1, cell] +
            soil_ice_available[1, cell] + soil_ice_free[1, cell]
        influx = compute_infiltration_slug(
            slug, top_storage, soil_wsats[1, cell] - soil_wpwps[1, cell], soil_infil,
        )
        soil_w_influx[1, cell] += influx
        srunoff = slug-influx
        soil_srunoff[cell] += slug - influx # surface runoff used for leaching
        incoming_volumetric_enthalpy = top_volumetric_enthalpy

        if transfer_heat
            rain_energy = influx * T(0.001) * rain_fraction * rain_volumetric_enthalpy
            melt_energy = influx * T(0.001) * melt_fraction * melt_volumetric_enthalpy
            rain_energy_input[cell] += rain_energy
            snowmelt_energy_input[cell] += melt_energy
        end

        for l in 1:soil_layers
            lrunoff = zero(T)
            layer_influx = influx
            influx = zero(T)

            if transfer_heat
                soil_perc_energy[l, cell] +=
                    layer_influx * T(0.001) * incoming_volumetric_enthalpy
            end

            # Nitrate percolated from the layer above enters this layer even
            # when water does not continue percolating out of it today.
            soil_NO3[l, cell] += NO3perc_ly
            NO3perc_ly = zero(T)

            soil_w[l, cell] += (soil_w_fw[l, cell] + layer_influx) / soil_whcs[l, cell]
            soil_w_fw[l, cell] = zero(T)

            layer_volumetric_enthalpy = compute_liquid_water_enthalpy(
                soil_temperature[l, cell], water_heat_capacity, ice_heat_capacity,
                volumetric_fusion_heat,
            )

            # Handle lateral runoff of water above saturation
            liquid_capacity = max(
                soil_layer_depth[l] - soil_freeze_depth[l, cell], zero(T),
            ) * (soil_wsat[l, cell] - soil_wpwp[l, cell])
            if (soil_w[l, cell] * soil_whcs[l, cell]) > liquid_capacity
                grunoff = (soil_w[l, cell] * soil_whcs[l, cell]) - liquid_capacity
                soil_w[l, cell] -= grunoff / soil_whcs[l, cell]
                soil_lrunoff[l, cell] += grunoff
                lrunoff += grunoff
            end

            # Additional saturation check
            if (soil_wpwps[l, cell] + soil_w[l, cell] * soil_whcs[l, cell] + soil_ice_available[l, cell] + soil_ice_free[l, cell]) > soil_wsats[l, cell]
                grunoff = (soil_wpwps[l, cell] + soil_w[l, cell] * soil_whcs[l, cell] + soil_ice_available[l, cell] + soil_ice_free[l, cell]) - soil_wsats[l, cell]
                soil_w[l, cell] -= grunoff / soil_whcs[l, cell]
                soil_lrunoff[l, cell] += grunoff
                lrunoff += grunoff
            end

            if transfer_heat && lrunoff > zero(T)
                runoff_energy = lrunoff * T(0.001) * layer_volumetric_enthalpy
                soil_perc_energy[l, cell] -= runoff_energy
                lateral_energy_output[cell] += runoff_energy
            end

            # Percolation from layer l to l+1 (or to outflux at bottom layer).
            if (soil_w[l, cell] - percthres) > (1.0f-5 / soil_whcs[l, cell])
                ice_fraction = soil_wsats[l, cell] > eps(T) ?
                               clamp(soil_ice[l, cell] / soil_wsats[l, cell], zero(T), one(T)) : zero(T)
                if l < soil_layers
                    saturation_factor = 1 - (soil_w[l+1, cell] * soil_whcs[l+1, cell] + soil_w_fw[l+1, cell] + soil_ice_available[l+1, cell] + soil_ice_free[l+1, cell]) / (soil_wsats[l+1, cell] - soil_wpwps[l+1, cell])
                else
                    saturation_factor = 1 - (soil_w[l, cell] * soil_whcs[l, cell] + soil_w_fw[l, cell] + soil_ice_available[l, cell] + soil_ice_free[l, cell]) / (soil_wsats[l, cell] - soil_wpwps[l, cell])
                end
                saturation = (soil_w[l, cell] * soil_whcs[l, cell] + soil_wpwps[l, cell] +
                    soil_ice_available[l, cell] + soil_ice_free[l, cell]) / soil_wsats[l, cell]
                perc = compute_percolation(
                    soil_w[l, cell], soil_whcs[l, cell], percthres, soil_Ks[l, cell],
                    saturation, soil_beta_soil[l, cell], ice_fraction, saturation_factor,
                )

                soil_w[l, cell] -= perc / soil_whcs[l, cell]

                if soil_w[l, cell] < 0
                    perc += soil_w[l, cell] * soil_whcs[l, cell]
                    soil_w[l, cell] = zero(T)
                end

                if l == soil_layers
                    soil_outflux_f[cell] += perc
                    soil_w_outflux[l, cell] += perc
                    if transfer_heat
                        drainage_energy = perc * T(0.001) * layer_volumetric_enthalpy
                        soil_perc_energy[l, cell] -= drainage_energy
                        bottom_energy_output[cell] += drainage_energy
                    end
                else
                    influx = perc
                    soil_w_influx[l+1, cell] += perc
                    soil_w_outflux[l, cell] += perc
                    if transfer_heat
                        transfer_energy = perc * T(0.001) * layer_volumetric_enthalpy
                        soil_perc_energy[l, cell] -= transfer_energy
                        # The next layer receives precisely the enthalpy carried
                        # out of this layer, not the upper-boundary mixture.
                        incoming_volumetric_enthalpy = layer_volumetric_enthalpy
                    end
                end

                w_mobile = perc + srunoff + lrunoff
                concNO3_mobile = compute_mobile_nitrate_concentration(
                    soil_NO3[l, cell], w_mobile, anion_excl, soil_wsats[l, cell],
                )
                # Surface runoff NO3 can be added here if a dedicated pool is introduced.
                srunoff = zero(T)
                NO3lat = zero(T)
                if l == 1
                    NO3lat = NPERCO * concNO3_mobile * lrunoff
                else
                    NO3lat = concNO3_mobile * lrunoff
                end
                NO3lat = min(NO3lat, soil_NO3[l, cell])
                soil_NO3[l, cell] -= NO3lat
                soil_n_leaching[cell] += NO3lat

                # nitrate percolating from this layer
                NO3perc_ly = concNO3_mobile * perc
                NO3perc_ly = min(NO3perc_ly, soil_NO3[l, cell])
                soil_NO3[l, cell] -= NO3perc_ly
                if l == soil_layers
                    soil_n_leaching[cell] += NO3perc_ly
                end
            end
        end
    end

    # Preserve LPJmL's MAXITER fallback. Do not leave rainfall in the mutable
    # infiltration input because it would be omitted from the day's water
    # balance and silently disappear at the next process update.
    if iter == 1000 && infil[cell] > zero(T)
        soil_srunoff[cell] += infil[cell]
        infil[cell] = zero(T)
    end

    for l in 1:soil_layers
        # Lateral runoff can draw from water already stored in a layer when
        # freezing reduces its liquid pore capacity, so it must be included in
        # the absolute layer-stock update as well as in the outgoing ledger.
        soil_perc[l, cell] = soil_w_influx[l, cell] -
                             soil_w_outflux[l, cell] -
                             soil_lrunoff[l, cell]
    end

    # Rainfall progressively settles the tilled topsoil back toward its
    # untilled bulk density. LPJmL applies this after the day's infiltration,
    # so the updated factor affects hydraulic properties on the next day.
    top_infiltration = soil_w_influx[1, cell]
    settling_fraction = compute_tillage_settling_fraction(
        top_infiltration, soil_sand_fraction[1, cell], soil_layer_depth[1],
    )
    tillage_density_factor[1, cell] += settling_fraction *
        (one(T) - tillage_density_factor[1, cell])

end
