"""
nitrogen_transform!(soil; air_temperature=nothing, wind_speed=nothing,
                    lpjmlparams=lpjmlparams)

Apply the LPJmL-style daily mineralization/immobilization, nitrification,
denitrification, and NH₃-volatilization sequence. Internal and boundary
fluxes are stored layer-wise in `soil.nitrogen` for diagnostics.
"""
function nitrogen_transform!(soil;
                             air_temperature = nothing,
                             wind_speed = nothing,
    lpjmlparams::LPJmLParams = lpjmlparams)
    mineralize_nitrify!(soil; lpjmlparams = lpjmlparams)
    post_crop_nitrogen_losses!(
        soil;
        air_temperature = air_temperature,
        wind_speed = wind_speed,
        lpjmlparams = lpjmlparams,
    )
    return nothing
end

"""
    mineralize_nitrify!(soil; lpjmlparams=lpjmlparams)

Release SOM and litter nitrogen, satisfy litter immobilization demand, and
nitrify the remaining ammonium. LPJmL performs this stage before the daily
plant processes, so the newly mineralized nitrogen is available for same-day
crop uptake.
"""
function mineralize_nitrify!(soil;
                             lpjmlparams::LPJmLParams = lpjmlparams,
                             shift_fast = soil_decomposition_input(soil).shift_fast,
                             shift_slow = soil_decomposition_input(soil).shift_slow)
    soil_layers = size(soil_nitrogen_prognostic(soil).nitrate, 1)
    launch_custom!(
        mineralize_immobilize_kernel!,
        soil_carbon_fluxes(soil).decomposed_litter,
        size(soil_carbon_fluxes(soil).decomposed_litter, 2),
        soil_nitrogen_fluxes(soil).decomposed_litter,
        soil_nitrogen_fluxes(soil).decomposed_fast,
        soil_nitrogen_fluxes(soil).decomposed_slow,
        shift_fast,
        shift_slow,
        soil_nitrogen_prognostic(soil).ammonium,
        soil_nitrogen_prognostic(soil).nitrate,
        soil_nitrogen_prognostic(soil).fast,
        soil_nitrogen_prognostic(soil).slow,
        soil_properties(soil).layer_depth,
        soil_nitrogen_fluxes(soil).mineralization,
        soil_nitrogen_fluxes(soil).immobilization,
        kernel_constant(
            soil_carbon_fluxes(soil).decomposed_litter,
            (; lpjmlparams, soil_layers),
        ),
    )

    launch_1D!(
        nitrify_kernel!,
        soil_properties(soil).ph,
        soil_nitrogen_prognostic(soil).ammonium,
        soil_nitrogen_prognostic(soil).nitrate,
        soil_water_auxiliary(soil).relative_content,
        soil_water_auxiliary(soil).holding_capacity_storage,
        soil_water_auxiliary(soil).wilting_storage,
        soil_water_prognostic(soil).wilting_ice_fraction,
        soil_water_auxiliary(soil).free_water,
        soil_water_auxiliary(soil).saturation_storage,
        soil_thermal_prognostic(soil).temperature,
        soil_nitrogen_fluxes(soil).nitrification,
        soil_nitrogen_fluxes(soil).n2o_nitrification,
        kernel_constant(
            soil_properties(soil).ph, (; lpjmlparams, soil_layers),
        ),
    )

    return nothing
end

"""
    post_crop_nitrogen_losses!(soil; air_temperature=nothing,
                               wind_speed=nothing, lpjmlparams=lpjmlparams)

Apply denitrification and ammonia volatilization after crop uptake, matching
LPJmL's daily stand ordering.
"""
function post_crop_nitrogen_losses!(soil;
                                    air_temperature = nothing,
                                    wind_speed = nothing,
                                    lpjmlparams::LPJmLParams = lpjmlparams)
    soil_layers = size(soil_nitrogen_prognostic(soil).nitrate, 1)

    launch_1D!(
        denitrify_kernel!,
        soil_properties(soil).ph,
        soil_carbon_prognostic(soil).fast,
        soil_carbon_prognostic(soil).slow,
        soil_water_auxiliary(soil).relative_content,
        soil_water_auxiliary(soil).holding_capacity_storage,
        soil_water_auxiliary(soil).wilting_storage,
        soil_water_prognostic(soil).wilting_ice_fraction,
        soil_water_auxiliary(soil).free_water,
        soil_water_auxiliary(soil).saturation_storage,
        soil_thermal_prognostic(soil).temperature,
        soil_nitrogen_prognostic(soil).nitrate,
        soil_nitrogen_fluxes(soil).denitrification,
        soil_nitrogen_fluxes(soil).n2o_denitrification,
        soil_nitrogen_fluxes(soil).n2_denitrification,
        kernel_constant(
            soil_properties(soil).ph, (; lpjmlparams, soil_layers),
        ),
    )

    volatilization_temperature = air_temperature === nothing ?
        vec(@view(soil_thermal_prognostic(soil).temperature[1, :])) : air_temperature
    volatilization_wind = if wind_speed === nothing
        fallback = soil_decomposition_workspace(soil).surface_scratch_1
        fill!(fallback, eltype(fallback)(lpjmlparams.volatil_wind))
        fallback
    else
        wind_speed
    end
    launch_1D!(
        volatilization_kernel!,
        soil_properties(soil).ph,
        soil_nitrogen_prognostic(soil).ammonium,
        volatilization_temperature,
        volatilization_wind,
        soil_properties(soil).layer_depth,
        soil_nitrogen_fluxes(soil).volatilization,
        kernel_constant(soil_properties(soil).ph, lpjmlparams),
    )
    return nothing
end

"""
    compute_water_filled_pore_space(relative_water, holding, wilting, ice_fraction,
                                    free_water, saturation)

Return liquid water-filled pore space in `[0, 1]`. This is shared by the
LPJmL nitrification and denitrification response formulations; ice is removed
from the wilting-water contribution before the ratio is formed.
"""
@inline function compute_water_filled_pore_space(relative_water::T,
                                                  holding::T,
                                                  wilting::T,
                                                  wilting_ice_fraction::T,
                                                  free_water::T,
                                                  saturation::T) where {T <: AbstractFloat}
    liquid_water = relative_water * holding +
        wilting * (one(T) - wilting_ice_fraction) + free_water
    return clamp(liquid_water / max(saturation, eps(T)), zero(T), one(T))
end

"""
    compute_nitrification_moisture_response(wfps, b, c, d, n, m, z)

Evaluate LPJmL's bounded two-branch water-filled-pore-space response. The
explicit positive-base guard avoids fractional powers of a negative value.
"""
@inline function compute_nitrification_moisture_response(wfps::T,
                                                          b::T,
                                                          c::T,
                                                          d::T,
                                                          n::T,
                                                          m::T,
                                                          z::T) where {T <: AbstractFloat}
    first_base = (wfps - b) / n
    second_base = (wfps - c) / m
    (first_base > zero(T) && second_base > zero(T)) || return zero(T)
    return max(zero(T), first_base^z * second_base^d)
end

"""
    compute_nitrification_temperature_response(temperature)

Gaussian optimum-temperature response used by LPJmL nitrification.
"""
@inline function compute_nitrification_temperature_response(temperature::T) where {T <: AbstractFloat}
    return exp(-(temperature - T(18.79))^2 / T(2 * 8.26 * 8.26))
end

"""
    compute_nitrification_ph_response(ph)

Smooth pH multiplier for nitrification. It intentionally remains unclamped to
match the source formulation and its calibrated parameter range.
"""
@inline function compute_nitrification_ph_response(ph::T) where {T <: AbstractFloat}
    return T(0.56) + atan(T(pi) * T(0.45) * (ph - T(5))) / T(pi)
end

"""
    compute_denitrification_temperature_response(temperature)

Piecewise LPJmL denitrification temperature response, including the low-
temperature baseline and high-temperature cutoff.
"""
@inline function compute_denitrification_temperature_response(temperature::T) where {T <: AbstractFloat}
    temperature > T(45.9) && return zero(T)
    temperature <= zero(T) && return T(0.0326)
    return max(zero(T), T(0.0326) + T(0.00351) * temperature^T(1.652) -
        (temperature / T(41.748))^T(7.19))
end

"""
    compute_denitrification_moisture_response(wfps)

Exponential anaerobic moisture limitation capped at one.
"""
@inline function compute_denitrification_moisture_response(wfps::T) where {T <: AbstractFloat}
    return min(one(T), T(6.664096e-10) * exp(T(20.92912) * wfps))
end

"""
    compute_ammonia_volatilization(ammonium, ph, temperature, wind, depth, length)

Return the LPJmL daily NH₃ loss from the upper ammonium pool. All quantities
are scalar, so the balance-preserving subtraction remains visible in the
calling kernel.
"""
@inline function compute_ammonia_volatilization(ammonium::T,
                                                 ph::T,
                                                 temperature::T,
                                                 wind::T,
                                                 depth::T,
                                                 volatilization_length::T) where {T <: AbstractFloat}
    available = max(zero(T), ammonium)
    kelvin = temperature + T(273.15)
    dissociation = T(10)^(T(0.05) - T(2788) / kelvin)
    aqueous_fraction = one(T) / (one(T) + T(10)^(-ph) / max(dissociation, eps(T)))
    aqueous_nh3 = aqueous_fraction * available / max(depth, eps(T)) * T(1000)
    henry = T(0.2138) / kelvin * T(10)^(T(6.123) - T(1825) / kelvin)
    transfer = T(0.000612) * max(zero(T), wind)^T(0.8) * kelvin^T(0.382) *
        volatilization_length^T(-0.2)
    return clamp(T(86400) * transfer * henry * aqueous_nh3, zero(T), available)
end

@kernel inbounds = true function mineralize_immobilize_kernel!(
    decomposed_litter_carbon::AbstractMatrix{T},
    decomposed_litter_nitrogen::AbstractMatrix{T},
    decomposed_fast_nitrogen::AbstractArray{M},
    decomposed_slow_nitrogen::AbstractArray{M},
    shift_fast::AbstractArray{M},
    shift_slow::AbstractArray{M},
    ammonium::AbstractArray{M},
    nitrate::AbstractArray{M},
    fast_nitrogen::AbstractArray{M},
    slow_nitrogen::AbstractArray{M},
    layer_depth::AbstractArray{T},
    mineralization::AbstractArray{M},
    immobilization::AbstractArray{M},
    fixed_parameters,
) where {T <: AbstractFloat, M <: AbstractFloat}
    cell = @index(Global)
    kernel_params = kernel_value(fixed_parameters)
    @unpack lpjmlparams, soil_layers = kernel_params
    @unpack atmfrac, fastfrac, soil_cn_ratio, immobilization_k = lpjmlparams

    litter_carbon = decomposed_litter_carbon[1, cell] +
        decomposed_litter_carbon[2, cell] + decomposed_litter_carbon[3, cell]
    litter_nitrogen = decomposed_litter_nitrogen[1, cell] +
        decomposed_litter_nitrogen[2, cell] + decomposed_litter_nitrogen[3, cell]
    carbon_nitrogen_deficit = litter_carbon / T(soil_cn_ratio) - litter_nitrogen

    for layer in 1:soil_layers
        mineralization[layer, cell] = zero(M)
        immobilization[layer, cell] = zero(M)

        # LPJmL keeps c_shift as a normalized vertical distribution. Apply the
        # fast/slow split and atmospheric fraction explicitly to each flux.
        litter_mineralization = max(
            zero(M),
            litter_nitrogen * T(atmfrac) *
            (T(fastfrac) * shift_fast[layer, cell] +
             (one(T) - T(fastfrac)) * shift_slow[layer, cell]),
        )
        som_mineralization = max(
            zero(M),
            decomposed_fast_nitrogen[layer, cell] +
            decomposed_slow_nitrogen[layer, cell],
        )
        gross_mineralization = litter_mineralization + som_mineralization
        ammonium[layer, cell] += gross_mineralization
        mineralization[layer, cell] = gross_mineralization

        if carbon_nitrogen_deficit > zero(T)
            available = max(zero(M), ammonium[layer, cell] + nitrate[layer, cell])
            if available > zero(M)
                concentration = available / max(layer_depth[layer], eps(T)) * T(1000)
                limitation = concentration / (T(immobilization_k) + concentration)

                fast_immobilization = max(
                    zero(M),
                    carbon_nitrogen_deficit * T(fastfrac) *
                    (one(T) - T(atmfrac)) * shift_fast[layer, cell] * limitation,
                )
                fast_immobilization = min(fast_immobilization, available)
                if fast_immobilization > zero(M)
                    ammonium_share = ammonium[layer, cell] / available
                    ammonium[layer, cell] -= fast_immobilization * ammonium_share
                    nitrate[layer, cell] -= fast_immobilization * (one(M) - ammonium_share)
                    fast_nitrogen[layer, cell] += fast_immobilization
                    immobilization[layer, cell] += fast_immobilization
                end

                available = max(zero(M), ammonium[layer, cell] + nitrate[layer, cell])
                if available > zero(M)
                    concentration = available / max(layer_depth[layer], eps(T)) * T(1000)
                    limitation = concentration / (T(immobilization_k) + concentration)
                    slow_immobilization = max(
                        zero(M),
                        carbon_nitrogen_deficit * (one(T) - T(fastfrac)) *
                        (one(T) - T(atmfrac)) * shift_slow[layer, cell] * limitation,
                    )
                    slow_immobilization = min(slow_immobilization, available)
                    if slow_immobilization > zero(M)
                        ammonium_share = ammonium[layer, cell] / available
                        ammonium[layer, cell] -= slow_immobilization * ammonium_share
                        nitrate[layer, cell] -= slow_immobilization * (one(M) - ammonium_share)
                        slow_nitrogen[layer, cell] += slow_immobilization
                        immobilization[layer, cell] += slow_immobilization
                    end
                end
            end
        end
        ammonium[layer, cell] = max(zero(M), ammonium[layer, cell])
        nitrate[layer, cell] = max(zero(M), nitrate[layer, cell])
    end
end

@kernel inbounds = true function nitrify_kernel!(
    soil_ph::AbstractArray{T},
    ammonium::AbstractArray{M},
    nitrate::AbstractArray{M},
    relative_water::AbstractArray{M},
    holding_storage::AbstractArray{M},
    wilting_storage::AbstractArray{M},
    wilting_ice_fraction::AbstractArray{M},
    free_water::AbstractArray{M},
    saturation_storage::AbstractArray{M},
    soil_temperature::AbstractArray{M},
    nitrification::AbstractArray{M},
    n2o_nitrification::AbstractArray{M},
    fixed_parameters,
) where {T <: AbstractFloat, M <: AbstractFloat}
    cell = @index(Global)
    kernel_params = kernel_value(fixed_parameters)
    @unpack lpjmlparams, soil_layers = kernel_params
    @unpack k_max, k_2, nitrification_a, nitrification_b,
            nitrification_c, nitrification_d = lpjmlparams

    for layer in 1:soil_layers
        nitrification[layer, cell] = zero(M)
        n2o_nitrification[layer, cell] = zero(M)
        water_filled_pore_space = compute_water_filled_pore_space(
            relative_water[layer, cell], holding_storage[layer, cell],
            wilting_storage[layer, cell], wilting_ice_fraction[layer, cell],
            free_water[layer, cell], saturation_storage[layer, cell],
        )
        n_nit = M(nitrification_a - nitrification_b)
        m_nit = M(nitrification_a - nitrification_c)
        z_nit = M(nitrification_d) * M(nitrification_b - nitrification_a) /
            M(nitrification_a - nitrification_c)
        moisture_factor = compute_nitrification_moisture_response(
            water_filled_pore_space, M(nitrification_b), M(nitrification_c),
            M(nitrification_d), n_nit, m_nit, z_nit,
        )
        temperature_factor = compute_nitrification_temperature_response(
            soil_temperature[layer, cell],
        )
        ph_factor = compute_nitrification_ph_response(M(soil_ph[cell]))
        gross_nitrification = clamp(
            T(k_max) * ammonium[layer, cell] * temperature_factor *
            moisture_factor * ph_factor,
            zero(M), ammonium[layer, cell],
        )
        n2o_loss = T(k_2) * gross_nitrification
        ammonium[layer, cell] -= gross_nitrification
        nitrate[layer, cell] += gross_nitrification - n2o_loss
        nitrification[layer, cell] = gross_nitrification
        n2o_nitrification[layer, cell] = n2o_loss
    end
end

@kernel inbounds = true function denitrify_kernel!(
    _soil_ph::AbstractArray{T},
    fast_carbon::AbstractArray{M},
    slow_carbon::AbstractArray{M},
    relative_water::AbstractArray{M},
    holding_storage::AbstractArray{M},
    wilting_storage::AbstractArray{M},
    wilting_ice_fraction::AbstractArray{M},
    free_water::AbstractArray{M},
    saturation_storage::AbstractArray{M},
    soil_temperature::AbstractArray{M},
    nitrate::AbstractArray{M},
    denitrification::AbstractArray{M},
    n2o_denitrification::AbstractArray{M},
    n2_denitrification::AbstractArray{M},
    fixed_parameters,
) where {T <: AbstractFloat, M <: AbstractFloat}
    cell = @index(Global)
    kernel_params = kernel_value(fixed_parameters)
    @unpack lpjmlparams, soil_layers = kernel_params
    @unpack CDN, n2o_denit_frac = lpjmlparams

    for layer in 1:soil_layers
        denitrification[layer, cell] = zero(M)
        n2o_denitrification[layer, cell] = zero(M)
        n2_denitrification[layer, cell] = zero(M)
        temperature = soil_temperature[layer, cell]
        organic_carbon = max(zero(M), fast_carbon[layer, cell] + slow_carbon[layer, cell])
        temperature_factor = compute_denitrification_temperature_response(temperature)
        water_filled_pore_space = compute_water_filled_pore_space(
            relative_water[layer, cell], holding_storage[layer, cell],
            wilting_storage[layer, cell], wilting_ice_fraction[layer, cell],
            free_water[layer, cell], saturation_storage[layer, cell],
        )
        moisture_factor = compute_denitrification_moisture_response(
            water_filled_pore_space,
        )
        carbon_factor = max(
            zero(M), one(M) - exp(-M(CDN) * temperature_factor * organic_carbon),
        )
        gross_denitrification = clamp(
            moisture_factor * carbon_factor * nitrate[layer, cell],
            zero(M), nitrate[layer, cell],
        )
        n2o_loss = M(n2o_denit_frac) * gross_denitrification
        n2_loss = gross_denitrification - n2o_loss
        nitrate[layer, cell] -= gross_denitrification
        denitrification[layer, cell] = gross_denitrification
        n2o_denitrification[layer, cell] = n2o_loss
        n2_denitrification[layer, cell] = n2_loss
    end
end

@kernel inbounds = true function volatilization_kernel!(
    soil_ph::AbstractArray{T},
    ammonium::AbstractArray{M},
    air_temperature::AbstractArray{M},
    wind_speed::AbstractArray{M},
    layer_depth::AbstractArray{T},
    volatilization::AbstractArray{M},
    fixed_parameters,
) where {T <: AbstractFloat, M <: AbstractFloat}
    cell = @index(Global)
    lpjmlparams = kernel_value(fixed_parameters)
    @unpack volatil_length = lpjmlparams
    flux = compute_ammonia_volatilization(
        ammonium[1, cell], M(soil_ph[cell]), air_temperature[cell], wind_speed[cell],
        M(layer_depth[1]), M(volatil_length),
    )
    ammonium[1, cell] -= flux
    volatilization[cell] = flux
end
