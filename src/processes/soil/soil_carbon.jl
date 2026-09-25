"""
soil_carbon!(crop, soil)

Update litter and soil carbon pools and heterotrophic respiration terms.
"""

"""Decompose existing litter and SOM carbon without routing new harvest residues."""
function soil_carbon_decomposition!(soil;
                                    lpjmlparams::LPJmLParams = lpjmlparams,
                                    soil_decomp_params::SoilDecompParams = soil_decomp_params)
    soil_decomp_response!(
        soil; lpjmlparams = lpjmlparams, soil_decomp_params = soil_decomp_params,
    )
    _soil_carbon_decomposition_from_response!(soil; lpjmlparams = lpjmlparams)
    return nothing
end

"""Apply carbon turnover using response arrays already stored on `soil`."""
function _soil_carbon_decomposition_from_response!(
    soil;
    lpjmlparams::LPJmLParams = lpjmlparams,
)
    launch_custom!(
        soil_carbon_decomposition_kernel!,
        soil_carbon_prognostic(soil).litter,
        size(soil_carbon_prognostic(soil).litter, 2),
        soil_carbon_auxiliary(soil).litter_response,
        soil_decomposition_auxiliary(soil).litter_response,
        soil_carbon_fluxes(soil).decomposed_litter,
        soil_carbon_prognostic(soil).fast,
        soil_carbon_prognostic(soil).slow,
        soil_decomposition_auxiliary(soil).response,
        soil_carbon_fluxes(soil).decomposed_fast,
        soil_carbon_fluxes(soil).decomposed_slow,
        soil_decomposition_input(soil).shift_fast,
        soil_decomposition_input(soil).shift_slow,
        soil_carbon_fluxes(soil).litter_to_fast,
        soil_carbon_fluxes(soil).litter_to_slow,
        soil_carbon_fluxes(soil).heterotrophic_respiration,
        kernel_constant(soil_carbon_prognostic(soil).litter, lpjmlparams),
        size(soil_carbon_prognostic(soil).fast, 1),
    )
    return nothing
end

"""
    soil_carbon!(crop, soil; ...)

Compatibility entry point for the former combined operation: decompose the
carbon pools, then route harvest-day residues without decomposing those new
residues until the following day.
"""
function soil_carbon!(crop,
                      soil;
                      lpjmlparams::LPJmLParams = lpjmlparams,
                      soil_decomp_params::SoilDecompParams = soil_decomp_params)
    soil_carbon_decomposition!(
        soil; lpjmlparams = lpjmlparams, soil_decomp_params = soil_decomp_params,
    )
    route_harvest_carbon_input!(soil, crop)
    return nothing
end

"""
    compute_first_order_decomposition(pool, rate, environmental_response)

Exact one-day first-order pool loss. `-expm1(-x)` is numerically stable for
small daily rates. Pool-specific non-negative guards remain explicit in the
calling kernel, preserving the original litter/SOM semantics.
"""
@inline function compute_first_order_decomposition(pool::T,
                                                   rate::T,
                                                   environmental_response::T) where {T <: AbstractFloat}
    return -expm1(-rate * environmental_response) * pool
end

"""
    compute_litter_to_som_routing(litter_flux, vertical_fraction,
                                  atmospheric_fraction, pool_fraction)

Return the non-respiratory share of decomposed litter delivered to one SOM
pool and one layer. `pool_fraction` selects the fast or slow partition.
"""
@inline function compute_litter_to_som_routing(litter_flux::T,
                                                vertical_fraction::T,
                                                atmospheric_fraction::T,
                                                pool_fraction::T) where {T <: AbstractFloat}
    return vertical_fraction * litter_flux * pool_fraction *
        (one(T) - atmospheric_fraction)
end

@kernel inbounds = true function soil_carbon_decomposition_kernel!(
    litter::AbstractMatrix{T},
    litter_rate::AbstractVector{T},
    litter_environment::AbstractMatrix{T},
    decomposed_litter::AbstractMatrix{T},
    fast::AbstractMatrix{T},
    slow::AbstractMatrix{T},
    soil_environment::AbstractMatrix{T},
    decomposed_fast::AbstractMatrix{T},
    decomposed_slow::AbstractMatrix{T},
    shift_fast::AbstractMatrix{T},
    shift_slow::AbstractMatrix{T},
    litter_to_fast::AbstractMatrix{T},
    litter_to_slow::AbstractMatrix{T},
    heterotrophic_respiration::AbstractVector{T},
    fixed_parameters,
    soil_layers::Integer,
) where {T <: AbstractFloat}
    cell = @index(Global)
    lpjmlparams = kernel_value(fixed_parameters)
    atmospheric_fraction = T(lpjmlparams.atmfrac)
    fast_fraction = T(lpjmlparams.fastfrac)
    fast_rate = T(lpjmlparams.k_soil10.fast)
    slow_rate = T(lpjmlparams.k_soil10.slow)
    litter_flux = zero(T)
    for pool in 1:3
        decomposition = compute_first_order_decomposition(
            litter[pool, cell], litter_rate[pool], litter_environment[pool, cell],
        )
        decomposed_litter[pool, cell] = decomposition
        litter[pool, cell] -= decomposition
        litter_flux += decomposition
    end

    respiration = litter_flux * atmospheric_fraction
    for layer in 1:soil_layers
        fast_decomposition = max(zero(T), compute_first_order_decomposition(
            fast[layer, cell], fast_rate, soil_environment[layer, cell],
        ))
        slow_decomposition = max(zero(T), compute_first_order_decomposition(
            slow[layer, cell], slow_rate, soil_environment[layer, cell],
        ))
        to_fast = compute_litter_to_som_routing(
            litter_flux, shift_fast[layer, cell], atmospheric_fraction, fast_fraction,
        )
        to_slow = compute_litter_to_som_routing(
            litter_flux, shift_slow[layer, cell], atmospheric_fraction,
            one(T) - fast_fraction,
        )
        decomposed_fast[layer, cell] = fast_decomposition
        decomposed_slow[layer, cell] = slow_decomposition
        litter_to_fast[layer, cell] = to_fast
        litter_to_slow[layer, cell] = to_slow
        fast[layer, cell] += to_fast - fast_decomposition
        slow[layer, cell] += to_slow - slow_decomposition
        respiration += fast_decomposition + slow_decomposition
    end
    heterotrophic_respiration[cell] = respiration
end
