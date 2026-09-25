"""
    soil_temperature!(soil, air_temperature, initial_temperature=air_temperature)

Advance the five-layer soil column by one day using an LPJmL-style enthalpy
formulation. Enthalpy below zero represents frozen soil, enthalpy between zero
and the volumetric fusion heat represents a 0 °C phase mixture, and enthalpy
above the fusion heat represents unfrozen soil. A fixed number of explicit
substeps keeps the nonlinear phase-change update predictable on CPU and GPU.

Snow and surface litter contribute serial thermal resistance. The kernel also
partitions each layer's conserved water stock into liquid water and ice and
records the daily surface-energy input and numerical energy residual.
"""
function soil_temperature!(soil,
                           air_temperature::AbstractArray{T},
                           initial_temperature::AbstractArray{T} = air_temperature;
                           thermalparams::SoilThermalParams{T} = soil_thermal_params,
                           snowparams::SnowParams{T} = snowparams) where {T <: AbstractFloat}
    launch_1D!(
        soil_temperature_kernel!,
        soil_thermal_prognostic(soil).initialized,
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
        soil_thermal_fluxes(soil).surface_energy_flux,
        soil_thermal_fluxes(soil).energy_residual,
        soil_thermal_fluxes(soil).untracked_water_energy_flux,
        soil_thermal_input(soil).diffusivity_0,
        soil_water_prognostic(soil).storage,
        soil_water_prognostic(soil).ice_storage,
        soil_water_auxiliary(soil).saturation_storage,
        soil_properties(soil).layer_depth,
        soil_snow_prognostic(soil).height,
        soil_surface_litter_prognostic(soil).depth,
        soil_surface_litter_prognostic(soil).water_storage,
        soil_surface_litter_prognostic(soil).temperature,
        soil_surface_litter_prognostic(soil).conductivity,
        air_temperature,
        initial_temperature,
        kernel_constant(
            soil_thermal_prognostic(soil).initialized,
            (; thermalparams, snowparams),
        ),
    )
    # Keep LPJmL's three liquid/ice reservoirs in a small second kernel. This
    # preserves the exact process order while avoiding a monolithic GPU kernel.
    partition_soil_water_ice!(soil)
    return nothing
end

@inline function enthalpy_temperature(enthalpy::T,
                                      heat_capacity_frozen::T,
                                      heat_capacity_unfrozen::T,
                                      latent_heat::T) where {T <: AbstractFloat}
    if enthalpy < zero(T)
        return enthalpy / max(heat_capacity_frozen, eps(T))
    elseif enthalpy > latent_heat
        return (enthalpy - latent_heat) /
               max(heat_capacity_unfrozen, eps(T))
    else
        return zero(T)
    end
end

@inline function enthalpy_frozen_fraction(enthalpy::T,
                                          latent_heat::T) where {T <: AbstractFloat}
    if latent_heat <= eps(T)
        return enthalpy <= zero(T) ? one(T) : zero(T)
    end
    return one(T) - clamp(enthalpy / latent_heat, zero(T), one(T))
end

@inline function temperature_phase_enthalpy(temperature::T,
                                            frozen_fraction::T,
                                            heat_capacity_frozen::T,
                                            heat_capacity_unfrozen::T,
                                            latent_heat::T) where {T <: AbstractFloat}
    if temperature < zero(T)
        return temperature * heat_capacity_frozen
    elseif temperature > zero(T)
        return temperature * heat_capacity_unfrozen + latent_heat
    else
        return (one(T) - clamp(frozen_fraction, zero(T), one(T))) * latent_heat
    end
end

@inline function harmonic_mean(a::T, b::T) where {T <: AbstractFloat}
    return T(2) * a * b / max(a + b, eps(T))
end

@inline function layer_thermal_properties(total_water::T,
                                          saturation_storage::T,
                                          layer_depth::T,
                                          diffusivity_0::T,
                                          params::SoilThermalParams{T}) where {T <: AbstractFloat}
    @unpack seconds_per_day, diffusivity_conversion, soil_heat_capacity,
            mineral_heat_capacity, water_heat_capacity, ice_heat_capacity,
            volumetric_fusion_heat, solid_conductivity, water_conductivity,
            ice_conductivity = params

    depth = max(layer_depth, eps(T))
    pore_fraction = clamp(saturation_storage / depth, zero(T), one(T))
    solid_storage = max(depth - saturation_storage, zero(T))
    water = max(total_water, zero(T))
    saturation = saturation_storage > eps(T) ?
                 clamp(water / saturation_storage, zero(T), one(T)) : zero(T)

    heat_capacity_frozen =
        (mineral_heat_capacity * solid_storage + ice_heat_capacity * water) / depth
    heat_capacity_unfrozen =
        (mineral_heat_capacity * solid_storage + water_heat_capacity * water) / depth
    latent_heat = water / depth * volumetric_fusion_heat

    dry_conductivity = max(
        diffusivity_0 * diffusivity_conversion / seconds_per_day *
        soil_heat_capacity,
        T(0.025),
    )
    saturated_frozen = T(10) ^ (
        log10(solid_conductivity) * (one(T) - pore_fraction) +
        log10(ice_conductivity) * pore_fraction
    )
    saturated_unfrozen = T(10) ^ (
        log10(solid_conductivity) * (one(T) - pore_fraction) +
        log10(water_conductivity) * pore_fraction
    )
    kersten_frozen = saturation
    kersten_unfrozen = saturation < T(0.1) ? zero(T) :
                       clamp(log10(saturation) + one(T), zero(T), one(T))
    conductivity_frozen = dry_conductivity +
        (saturated_frozen - dry_conductivity) * kersten_frozen
    conductivity_unfrozen = dry_conductivity +
        (saturated_unfrozen - dry_conductivity) * kersten_unfrozen

    return heat_capacity_frozen, heat_capacity_unfrozen, latent_heat,
           max(conductivity_frozen, T(0.025)),
           max(conductivity_unfrozen, T(0.025))
end

"""
    compute_litter_thermal_conductivity(water_storage, depth, temperature, params)

Evaluate the dry-to-saturated litter conductivity interpolation used by the
surface thermal resistance. This is kept separate from the five-layer heat
solver because it is a local constitutive relation, not a column update.
"""
@inline function compute_litter_thermal_conductivity(
    water_storage::T,
    depth::T,
    temperature::T,
    porosity::T,
    dry_conductivity::T,
    saturated_unfrozen::T,
    saturated_frozen::T,
) where {T <: AbstractFloat}
    depth > eps(T) || return dry_conductivity
    saturation = water_storage > eps(T) ?
        min((water_storage / T(1000)) / (porosity * depth), one(T)) : zero(T)
    saturation <= eps(T) && return dry_conductivity
    if temperature < zero(T)
        return dry_conductivity + (saturated_frozen - dry_conductivity) * saturation
    end
    kersten = saturation < T(0.1) ? zero(T) : log10(saturation) + one(T)
    return dry_conductivity + (saturated_unfrozen - dry_conductivity) * kersten
end

@kernel inbounds = true function soil_temperature_kernel!(
    initialized::AbstractArray{Bool},
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
    surface_energy_flux::AbstractArray{T},
    energy_residual::AbstractArray{T},
    untracked_water_energy_flux::AbstractArray{T},
    diffusivity_0::AbstractArray{T},
    liquid_water::AbstractArray{T},
    ice_water::AbstractArray{T},
    saturation_storage::AbstractArray{T},
    layer_depth_mm::AbstractArray{T},
    snow_height::AbstractArray{T},
    litter_depth::AbstractArray{T},
    litter_water_storage::AbstractArray{T},
    litter_temperature::AbstractArray{T},
    litter_conductivity::AbstractArray{T},
    air_temperature::AbstractArray{T},
    initial_temperature::AbstractArray{T},
    fixed_parameters,
) where {T <: AbstractFloat}
    cell = @index(Global)
    parameters = kernel_value(fixed_parameters)
    thermalparams = parameters.thermalparams
    snowparams = parameters.snowparams
    snow_conductivity = T(snowparams.lambda_snow)
    @unpack seconds_per_day, phase_change_substeps, litter_porosity,
            litter_conductivity_dry, litter_conductivity_saturated_unfrozen,
            litter_conductivity_saturated_frozen = thermalparams

    dz1 = layer_depth_mm[1] * T(0.001)
    dz2 = layer_depth_mm[2] * T(0.001)
    dz3 = layer_depth_mm[3] * T(0.001)
    dz4 = layer_depth_mm[4] * T(0.001)
    dz5 = layer_depth_mm[5] * T(0.001)
    previous_e1 = enthalpy[1, cell]
    previous_e2 = enthalpy[2, cell]
    previous_e3 = enthalpy[3, cell]
    previous_e4 = enthalpy[4, cell]
    previous_e5 = enthalpy[5, cell]

    total1 = max(liquid_water[1, cell] + ice_water[1, cell], zero(T))
    total2 = max(liquid_water[2, cell] + ice_water[2, cell], zero(T))
    total3 = max(liquid_water[3, cell] + ice_water[3, cell], zero(T))
    total4 = max(liquid_water[4, cell] + ice_water[4, cell], zero(T))
    total5 = max(liquid_water[5, cell] + ice_water[5, cell], zero(T))

    cf1, cu1, lh1, kf1, ku1 = layer_thermal_properties(total1, saturation_storage[1, cell], layer_depth_mm[1], diffusivity_0[cell], thermalparams)
    cf2, cu2, lh2, kf2, ku2 = layer_thermal_properties(total2, saturation_storage[2, cell], layer_depth_mm[2], diffusivity_0[cell], thermalparams)
    cf3, cu3, lh3, kf3, ku3 = layer_thermal_properties(total3, saturation_storage[3, cell], layer_depth_mm[3], diffusivity_0[cell], thermalparams)
    cf4, cu4, lh4, kf4, ku4 = layer_thermal_properties(total4, saturation_storage[4, cell], layer_depth_mm[4], diffusivity_0[cell], thermalparams)
    cf5, cu5, lh5, kf5, ku5 = layer_thermal_properties(total5, saturation_storage[5, cell], layer_depth_mm[5], diffusivity_0[cell], thermalparams)

    if initialized[cell]
        oldt1 = enthalpy_temperature(enthalpy[1, cell], heat_capacity_frozen[1, cell], heat_capacity_unfrozen[1, cell], latent_heat[1, cell])
        oldt2 = enthalpy_temperature(enthalpy[2, cell], heat_capacity_frozen[2, cell], heat_capacity_unfrozen[2, cell], latent_heat[2, cell])
        oldt3 = enthalpy_temperature(enthalpy[3, cell], heat_capacity_frozen[3, cell], heat_capacity_unfrozen[3, cell], latent_heat[3, cell])
        oldt4 = enthalpy_temperature(enthalpy[4, cell], heat_capacity_frozen[4, cell], heat_capacity_unfrozen[4, cell], latent_heat[4, cell])
        oldt5 = enthalpy_temperature(enthalpy[5, cell], heat_capacity_frozen[5, cell], heat_capacity_unfrozen[5, cell], latent_heat[5, cell])
        oldf1 = enthalpy_frozen_fraction(enthalpy[1, cell], latent_heat[1, cell])
        oldf2 = enthalpy_frozen_fraction(enthalpy[2, cell], latent_heat[2, cell])
        oldf3 = enthalpy_frozen_fraction(enthalpy[3, cell], latent_heat[3, cell])
        oldf4 = enthalpy_frozen_fraction(enthalpy[4, cell], latent_heat[4, cell])
        oldf5 = enthalpy_frozen_fraction(enthalpy[5, cell], latent_heat[5, cell])
        e1 = temperature_phase_enthalpy(oldt1, oldf1, cf1, cu1, lh1)
        e2 = temperature_phase_enthalpy(oldt2, oldf2, cf2, cu2, lh2)
        e3 = temperature_phase_enthalpy(oldt3, oldf3, cf3, cu3, lh3)
        e4 = temperature_phase_enthalpy(oldt4, oldf4, cf4, cu4, lh4)
        e5 = temperature_phase_enthalpy(oldt5, oldf5, cf5, cu5, lh5)
    else
        inittemp = initial_temperature[cell]
        initf1 = total1 > eps(T) ? clamp(ice_water[1, cell] / total1, zero(T), one(T)) : zero(T)
        initf2 = total2 > eps(T) ? clamp(ice_water[2, cell] / total2, zero(T), one(T)) : zero(T)
        initf3 = total3 > eps(T) ? clamp(ice_water[3, cell] / total3, zero(T), one(T)) : zero(T)
        initf4 = total4 > eps(T) ? clamp(ice_water[4, cell] / total4, zero(T), one(T)) : zero(T)
        initf5 = total5 > eps(T) ? clamp(ice_water[5, cell] / total5, zero(T), one(T)) : zero(T)
        e1 = temperature_phase_enthalpy(inittemp, initf1, cf1, cu1, lh1)
        e2 = temperature_phase_enthalpy(inittemp, initf2, cf2, cu2, lh2)
        e3 = temperature_phase_enthalpy(inittemp, initf3, cf3, cu3, lh3)
        e4 = temperature_phase_enthalpy(inittemp, initf4, cf4, cu4, lh4)
        e5 = temperature_phase_enthalpy(inittemp, initf5, cf5, cu5, lh5)
    end

    rebased_column_energy = e1 * dz1 + e2 * dz2 + e3 * dz3 + e4 * dz4 + e5 * dz5
    # LPJmL's apply_enth_of_untracked_mass_shifts: water changes not already
    # represented by perc_energy (primarily yesterday's evaporation and
    # transpiration) are assigned the enthalpy of water in the same layer, so
    # the mass change itself does not create a temperature jump.
    # Sum layer-wise changes instead of subtracting two O(1e8) column totals.
    # The expressions are algebraically equivalent, but this avoids losing the
    # small boundary flux to Float32 cancellation differently on CPU and GPU.
    untracked_water_energy_flux[cell] = initialized[cell] ?
        (e1 - previous_e1) * dz1 + (e2 - previous_e2) * dz2 +
        (e3 - previous_e3) * dz3 + (e4 - previous_e4) * dz4 +
        (e5 - previous_e5) * dz5 : zero(T)

    heat_capacity_frozen[1, cell] = cf1; heat_capacity_frozen[2, cell] = cf2
    heat_capacity_frozen[3, cell] = cf3; heat_capacity_frozen[4, cell] = cf4
    heat_capacity_frozen[5, cell] = cf5
    heat_capacity_unfrozen[1, cell] = cu1; heat_capacity_unfrozen[2, cell] = cu2
    heat_capacity_unfrozen[3, cell] = cu3; heat_capacity_unfrozen[4, cell] = cu4
    heat_capacity_unfrozen[5, cell] = cu5
    latent_heat[1, cell] = lh1; latent_heat[2, cell] = lh2
    latent_heat[3, cell] = lh3; latent_heat[4, cell] = lh4; latent_heat[5, cell] = lh5
    conductivity_frozen[1, cell] = kf1; conductivity_frozen[2, cell] = kf2
    conductivity_frozen[3, cell] = kf3; conductivity_frozen[4, cell] = kf4
    conductivity_frozen[5, cell] = kf5
    conductivity_unfrozen[1, cell] = ku1; conductivity_unfrozen[2, cell] = ku2
    conductivity_unfrozen[3, cell] = ku3; conductivity_unfrozen[4, cell] = ku4
    conductivity_unfrozen[5, cell] = ku5

    depth_litter = max(litter_depth[cell], zero(T))
    old_litter_temperature = initialized[cell] ? litter_temperature[cell] : initial_temperature[cell]
    conductivity_litter = compute_litter_thermal_conductivity(
        litter_water_storage[cell], depth_litter, old_litter_temperature,
        litter_porosity, litter_conductivity_dry,
        litter_conductivity_saturated_unfrozen, litter_conductivity_saturated_frozen,
    )
    litter_conductivity[cell] = conductivity_litter

    energy_before = rebased_column_energy
    surface_energy = zero(T)
    dt = seconds_per_day / T(phase_change_substeps)
    last_surface_resistance = one(T)
    for _ in 1:phase_change_substeps
        t1 = enthalpy_temperature(e1, cf1, cu1, lh1)
        t2 = enthalpy_temperature(e2, cf2, cu2, lh2)
        t3 = enthalpy_temperature(e3, cf3, cu3, lh3)
        t4 = enthalpy_temperature(e4, cf4, cu4, lh4)
        t5 = enthalpy_temperature(e5, cf5, cu5, lh5)
        ff1 = enthalpy_frozen_fraction(e1, lh1)
        ff2 = enthalpy_frozen_fraction(e2, lh2)
        ff3 = enthalpy_frozen_fraction(e3, lh3)
        ff4 = enthalpy_frozen_fraction(e4, lh4)
        ff5 = enthalpy_frozen_fraction(e5, lh5)
        k1 = ff1 * kf1 + (one(T) - ff1) * ku1
        k2 = ff2 * kf2 + (one(T) - ff2) * ku2
        k3 = ff3 * kf3 + (one(T) - ff3) * ku3
        k4 = ff4 * kf4 + (one(T) - ff4) * ku4
        k5 = ff5 * kf5 + (one(T) - ff5) * ku5

        surface_resistance = dz1 / (T(2) * k1) +
            depth_litter / conductivity_litter +
            max(snow_height[cell], zero(T)) / snow_conductivity
        qsurface = (air_temperature[cell] - t1) / max(surface_resistance, eps(T))
        q12 = harmonic_mean(k1, k2) * (t1 - t2) / ((dz1 + dz2) / T(2))
        q23 = harmonic_mean(k2, k3) * (t2 - t3) / ((dz2 + dz3) / T(2))
        q34 = harmonic_mean(k3, k4) * (t3 - t4) / ((dz3 + dz4) / T(2))
        q45 = harmonic_mean(k4, k5) * (t4 - t5) / ((dz4 + dz5) / T(2))
        e1 += dt * (qsurface - q12) / dz1
        e2 += dt * (q12 - q23) / dz2
        e3 += dt * (q23 - q34) / dz3
        e4 += dt * (q34 - q45) / dz4
        e5 += dt * q45 / dz5
        surface_energy += qsurface * dt
        last_surface_resistance = surface_resistance
    end

    enthalpy[1, cell] = e1; enthalpy[2, cell] = e2; enthalpy[3, cell] = e3
    enthalpy[4, cell] = e4; enthalpy[5, cell] = e5
    newt1 = enthalpy_temperature(e1, cf1, cu1, lh1)
    newt2 = enthalpy_temperature(e2, cf2, cu2, lh2)
    newt3 = enthalpy_temperature(e3, cf3, cu3, lh3)
    newt4 = enthalpy_temperature(e4, cf4, cu4, lh4)
    newt5 = enthalpy_temperature(e5, cf5, cu5, lh5)
    ff1 = enthalpy_frozen_fraction(e1, lh1); ff2 = enthalpy_frozen_fraction(e2, lh2)
    ff3 = enthalpy_frozen_fraction(e3, lh3); ff4 = enthalpy_frozen_fraction(e4, lh4)
    ff5 = enthalpy_frozen_fraction(e5, lh5)
    temperature[1, cell] = newt1; temperature[2, cell] = newt2
    temperature[3, cell] = newt3; temperature[4, cell] = newt4
    temperature[5, cell] = newt5
    frozen_fraction[1, cell] = ff1; frozen_fraction[2, cell] = ff2
    frozen_fraction[3, cell] = ff3; frozen_fraction[4, cell] = ff4
    frozen_fraction[5, cell] = ff5
    freeze_depth[1, cell] = ff1 * layer_depth_mm[1]
    freeze_depth[2, cell] = ff2 * layer_depth_mm[2]
    freeze_depth[3, cell] = ff3 * layer_depth_mm[3]
    freeze_depth[4, cell] = ff4 * layer_depth_mm[4]
    freeze_depth[5, cell] = ff5 * layer_depth_mm[5]

    ice_water[1, cell] = ff1 * total1; liquid_water[1, cell] = total1 - ice_water[1, cell]
    ice_water[2, cell] = ff2 * total2; liquid_water[2, cell] = total2 - ice_water[2, cell]
    ice_water[3, cell] = ff3 * total3; liquid_water[3, cell] = total3 - ice_water[3, cell]
    ice_water[4, cell] = ff4 * total4; liquid_water[4, cell] = total4 - ice_water[4, cell]
    ice_water[5, cell] = ff5 * total5; liquid_water[5, cell] = total5 - ice_water[5, cell]
    water_reference[1, cell] = total1; water_reference[2, cell] = total2
    water_reference[3, cell] = total3; water_reference[4, cell] = total4
    water_reference[5, cell] = total5

    energy_after = e1 * dz1 + e2 * dz2 + e3 * dz3 + e4 * dz4 + e5 * dz5
    surface_energy_flux[cell] = surface_energy
    energy_residual[cell] = energy_after - energy_before - surface_energy
    if depth_litter > eps(T)
        final_surface_flux = (air_temperature[cell] - newt1) /
                             max(last_surface_resistance, eps(T))
        litter_temperature[cell] = air_temperature[cell] - final_surface_flux *
            (max(snow_height[cell], zero(T)) / snow_conductivity +
             depth_litter / (T(2) * conductivity_litter))
    else
        litter_temperature[cell] = newt1
    end
    initialized[cell] = true
end
