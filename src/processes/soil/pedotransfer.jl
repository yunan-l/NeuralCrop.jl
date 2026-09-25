"""
pedotransfer!(soil; lpjmlparams=lpjmlparams)

Derive soil hydraulic properties from texture and depth parameterizations.
"""

"""Compute hydraulic properties in one backend-neutral layer/cell kernel."""
function pedotransfer!(soil;
                       lpjmlparams::LPJmLParams = lpjmlparams)
    launch_2D!(
        pedotransfer_kernel!,
        soil_water_auxiliary(soil).wilting_fraction,
        soil_properties(soil).sand_fraction,
        soil_properties(soil).clay_fraction,
        soil_properties(soil).layer_depth,
        soil_carbon_prognostic(soil).fast,
        soil_carbon_prognostic(soil).slow,
        soil_water_auxiliary(soil).wilting_storage,
        soil_water_auxiliary(soil).field_capacity,
        soil_water_prognostic(soil).saturation_fraction,
        soil_water_auxiliary(soil).saturation_storage,
        soil_water_auxiliary(soil).beta,
        soil_water_auxiliary(soil).holding_capacity_fraction,
        soil_water_auxiliary(soil).holding_capacity_storage,
        soil_water_auxiliary(soil).saturated_conductivity,
        soil_management_prognostic(soil).tillage_density_factor,
        kernel_constant(soil_water_auxiliary(soil).wilting_fraction, lpjmlparams),
    )
    partition_soil_water_ice!(soil)
    return nothing
end

function pedotransfer!(state::ModelState;
                       lpjmlparams::LPJmLParams = lpjmlparams)
    water_state = state.prognostic.soil.water
    water_auxiliary = state.auxiliary.soil.water
    soil_inputs = state.inputs.soil
    soil_carbon = state.prognostic.soil.carbon
    launch_2D!(
        pedotransfer_kernel!, water_auxiliary.wilting_fraction,
        soil_inputs.properties.sand_fraction,
        soil_inputs.properties.clay_fraction,
        soil_inputs.properties.layer_depth, soil_carbon.fast, soil_carbon.slow,
        water_auxiliary.wilting_storage, water_auxiliary.field_capacity,
        water_state.saturation_fraction, water_auxiliary.saturation_storage,
        water_auxiliary.beta, water_auxiliary.holding_capacity_fraction,
        water_auxiliary.holding_capacity_storage,
        water_auxiliary.saturated_conductivity,
        state.prognostic.soil.management.tillage_density_factor,
        kernel_constant(water_auxiliary.wilting_fraction, lpjmlparams),
    )
    partition_soil_water_ice!(state)
    return nothing
end

"""
    compute_soil_organic_matter(fast_carbon, slow_carbon, saturation, density, depth)

Convert layer carbon stocks to the bounded organic-matter percentage used by
the LPJmL/Saxton--Rawls pedotransfer relationships. This is a constitutive
relation only; it does not update any hydrological state.
"""
@inline function compute_soil_organic_matter(fast_carbon::T,
                                              slow_carbon::T,
                                              saturation::T,
                                              mineral_density::T,
                                              depth::T) where {T <: AbstractFloat}
    return clamp(
        T(2) * ((fast_carbon + slow_carbon) /
        ((one(T) - saturation) * mineral_density * depth)) * T(100),
        zero(T), T(8),
    )
end

"""
    compute_hydraulic_properties(sand, clay, organic_matter,
                                 tillage_density_factor, is_top_layer)

Evaluate the scalar Saxton--Rawls hydraulic parameterization for one soil
layer. The returned tuple is `(wilting, field, saturation, beta, holding,
conductivity)`; storage conversions remain in the kernel so the layer's array
updates and units stay explicit.
"""
@inline function compute_hydraulic_properties(sand::T,
                                              clay::T,
                                              organic_matter::T,
                                              tillage_density_factor::T,
                                              is_top_layer::Bool) where {T <: AbstractFloat}
    wpwpt = -T(0.024) * sand + T(0.487) * clay + T(0.006) * organic_matter +
        T(0.005) * sand * organic_matter - T(0.013) * clay * organic_matter +
        T(0.068) * sand * clay + T(0.031)
    wilting = wpwpt + (T(0.14) * wpwpt - T(0.02))

    ws33t = T(0.278) * sand + T(0.034) * clay + T(0.022) * organic_matter -
        T(0.018) * sand * organic_matter - T(0.027) * clay * organic_matter -
        T(0.584) * sand * clay + T(0.078)
    ws33 = ws33t + (T(0.636) * ws33t - T(0.107))

    wfct = -T(0.251) * sand + T(0.195) * clay + T(0.011) * organic_matter +
        T(0.006) * sand * organic_matter - T(0.027) * clay * organic_matter +
        T(0.452) * sand * clay + T(0.299)
    field = wfct + ((T(1.283) * wfct)^2 - T(0.374) * wfct - T(0.015))
    base_saturation = field + ws33 - T(0.097) * sand + T(0.043)
    saturation = is_top_layer ?
        one(T) - (one(T) - base_saturation) * tillage_density_factor :
        base_saturation
    is_top_layer && (field -= T(0.2) * (base_saturation - saturation))
    field = saturation - field < T(0.05) ? saturation - T(0.05) : field

    beta = -T(2.655) / log10(field / saturation)
    holding = field - wilting
    lambda = (log(field) - log(wilting)) / (log(T(1500)) - log(T(33)))
    conductivity = T(1930) * (saturation - field)^(T(3) - lambda)
    return wilting, field, saturation, beta, holding, conductivity
end

@kernel inbounds = true function pedotransfer_kernel!(
    wilting_fraction::AbstractMatrix{T},
    sand_fraction::AbstractMatrix{T},
    clay_fraction::AbstractMatrix{T},
    layer_depth::AbstractVector{T},
    fast_carbon::AbstractMatrix{T},
    slow_carbon::AbstractMatrix{T},
    wilting_storage::AbstractMatrix{T},
    field_capacity::AbstractMatrix{T},
    saturation_fraction::AbstractMatrix{T},
    saturation_storage::AbstractMatrix{T},
    beta::AbstractMatrix{T},
    holding_capacity_fraction::AbstractMatrix{T},
    holding_capacity_storage::AbstractMatrix{T},
    saturated_conductivity::AbstractMatrix{T},
    tillage_density_factor::AbstractMatrix{T},
    fixed_parameters,
) where {T <: AbstractFloat}
    layer, cell = @index(Global, NTuple)
    lpjmlparams = kernel_value(fixed_parameters)

    sand = sand_fraction[1, cell]
    clay = clay_fraction[1, cell]
    depth = layer_depth[layer]
    previous_saturation = saturation_fraction[layer, cell]
    organic_matter = compute_soil_organic_matter(
        fast_carbon[layer, cell], slow_carbon[layer, cell], previous_saturation,
        T(lpjmlparams.MINERALDENS), depth,
    )
    wilting, field, saturation, layer_beta, holding, conductivity =
        compute_hydraulic_properties(
            sand, clay, organic_matter, tillage_density_factor[1, cell], layer == 1,
        )

    wilting_fraction[layer, cell] = wilting
    wilting_storage[layer, cell] = wilting * depth
    field_capacity[layer, cell] = field
    saturation_fraction[layer, cell] = saturation
    saturation_storage[layer, cell] = saturation * depth
    beta[layer, cell] = layer_beta
    holding_capacity_fraction[layer, cell] = holding
    holding_capacity_storage[layer, cell] = holding * depth
    saturated_conductivity[layer, cell] = conductivity
end
