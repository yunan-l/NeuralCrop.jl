"""
soil_nitrogen!(crop, soil)

Update litter and soil nitrogen pools and crop-available mineral nitrogen.
"""

"""Decompose existing litter and SOM nitrogen without mineral transformations."""
function soil_nitrogen_decomposition!(soil;
                                      lpjmlparams::LPJmLParams = lpjmlparams,
                                      soil_decomp_params::SoilDecompParams = soil_decomp_params,
                                      litter_rate = soil_nitrogen_auxiliary(soil).litter_response,
                                      shift_fast = soil_decomposition_input(soil).shift_fast,
                                      shift_slow = soil_decomposition_input(soil).shift_slow)
    launch_custom!(
        soil_nitrogen_decomposition_kernel!,
        soil_nitrogen_prognostic(soil).litter,
        size(soil_nitrogen_prognostic(soil).litter, 2),
        litter_rate,
        soil_decomposition_auxiliary(soil).litter_response,
        soil_nitrogen_fluxes(soil).decomposed_litter,
        soil_nitrogen_prognostic(soil).fast,
        soil_nitrogen_prognostic(soil).slow,
        soil_decomposition_auxiliary(soil).response,
        soil_nitrogen_fluxes(soil).decomposed_fast,
        soil_nitrogen_fluxes(soil).decomposed_slow,
        shift_fast,
        shift_slow,
        soil_nitrogen_fluxes(soil).litter_to_fast,
        soil_nitrogen_fluxes(soil).litter_to_slow,
        kernel_constant(soil_nitrogen_prognostic(soil).litter, lpjmlparams),
        size(soil_nitrogen_prognostic(soil).fast, 1),
    )
    return nothing
end

"""
    soil_cn_decomposition!(soil; ...)

Execute LPJmL's coupled pre-crop soil stage: compute one shared environmental
response, decompose C and N with identical pool-specific decay fractions,
then mineralize/immobilize and nitrify mineral nitrogen.
"""
function soil_cn_decomposition!(soil;
                                lpjmlparams::LPJmLParams = lpjmlparams,
                                soil_decomp_params::SoilDecompParams = soil_decomp_params)
    soil_carbon_decomposition!(
        soil; lpjmlparams = lpjmlparams, soil_decomp_params = soil_decomp_params,
    )
    _soil_nitrogen_from_decomposition_response!(soil; lpjmlparams = lpjmlparams)
    return nothing
end

"""Complete coupled C--N turnover after a caller has supplied response arrays."""
function _soil_nitrogen_from_decomposition_response!(
    soil;
    lpjmlparams::LPJmLParams = lpjmlparams,
)
    soil_nitrogen_decomposition!(
        soil;
        lpjmlparams = lpjmlparams,
        litter_rate = soil_carbon_auxiliary(soil).litter_response,
        shift_fast = soil_decomposition_input(soil).shift_fast,
        shift_slow = soil_decomposition_input(soil).shift_slow,
    )
    mineralize_nitrify!(
        soil;
        lpjmlparams = lpjmlparams,
        shift_fast = soil_decomposition_input(soil).shift_fast,
        shift_slow = soil_decomposition_input(soil).shift_slow,
    )
    return nothing
end


"""Apply coupled C--N turnover using environmental responses already on `soil`."""
function _soil_cn_decomposition_from_response!(
    soil;
    lpjmlparams::LPJmLParams = lpjmlparams,
)
    _soil_carbon_decomposition_from_response!(soil; lpjmlparams = lpjmlparams)
    _soil_nitrogen_from_decomposition_response!(soil; lpjmlparams = lpjmlparams)
    return nothing
end

"""Route new harvest-day carbon and nitrogen residues together."""
function route_harvest_residues!(soil, crop)
    route_harvest_carbon_input!(soil, crop)
    route_harvest_nitrogen_input!(soil, crop)
    return nothing
end

"""
    soil_nitrogen!(crop, soil; ...)

Compatibility entry point for the former combined operation.
"""
function soil_nitrogen!(crop,
                        soil;
                        air_temperature = nothing,
                        wind_speed = nothing,
                        lpjmlparams::LPJmLParams = lpjmlparams,
                        soil_decomp_params::SoilDecompParams = soil_decomp_params)
    soil_nitrogen_decomposition!(
        soil; lpjmlparams = lpjmlparams, soil_decomp_params = soil_decomp_params,
    )
    route_harvest_nitrogen_input!(soil, crop)
    nitrogen_transform!(
        soil;
        air_temperature = air_temperature,
        wind_speed = wind_speed,
        lpjmlparams = lpjmlparams,
    )
    return nothing
end

@kernel inbounds = true function soil_nitrogen_decomposition_kernel!(
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
    end
end
