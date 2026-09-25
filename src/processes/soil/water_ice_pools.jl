"""
    lpjml_water_ice_partition(total_water, total_ice, wilting_storage,
                              holding_storage)

Partition a layer's conserved water and ice into LPJmL's three hydrological
reservoirs: water below the permanent wilting point, plant-available water,
and gravitational/free water. Ice is assigned to the first two reservoirs
with a common fraction; free water freezes only after both are fully frozen,
matching `freezefrac2soil()`.
"""
@inline function lpjml_water_ice_partition(total_water::T,
                                           total_ice::T,
                                           wilting_storage::T,
                                           holding_storage::T) where {T <: AbstractFloat}
    total = max(total_water, zero(T))
    ice = clamp(total_ice, zero(T), total)
    unavailable = min(total, max(wilting_storage, zero(T)))
    available = min(
        max(total - unavailable, zero(T)),
        max(holding_storage, zero(T)),
    )
    free = max(total - unavailable - available, zero(T))

    bound = unavailable + available
    bound_ice_fraction = bound > eps(T) ?
        min(ice / bound, one(T)) : zero(T)
    remaining_ice = max(ice - bound, zero(T))
    free_ice_fraction = free > eps(T) ?
        clamp(remaining_ice / free, zero(T), one(T)) : zero(T)

    unavailable_ice = bound_ice_fraction * unavailable
    available_ice = bound_ice_fraction * available
    free_ice = free_ice_fraction * free
    pwp_ice_fraction = wilting_storage > eps(T) ?
        clamp(unavailable_ice / wilting_storage, zero(T), one(T)) : zero(T)
    relative_water = holding_storage > eps(T) ?
        clamp((available - available_ice) / holding_storage, zero(T), one(T)) : zero(T)
    free_water = free - free_ice

    return pwp_ice_fraction, available_ice, free_ice,
           relative_water, free_water
end

"""Rebuild LPJmL liquid and ice reservoirs from conserved layer totals."""
function partition_soil_water_ice!(soil)
    launch_custom!(
        partition_soil_water_ice_kernel!,
        soil_water_prognostic(soil).storage,
        size(soil_water_prognostic(soil).storage, 2),
        soil_water_prognostic(soil).ice_storage,
        soil_water_prognostic(soil).wilting_ice_fraction,
        soil_water_prognostic(soil).available_ice_storage,
        soil_water_prognostic(soil).free_ice_storage,
        soil_water_auxiliary(soil).relative_content,
        soil_water_auxiliary(soil).free_water,
        soil_water_auxiliary(soil).wilting_storage,
        soil_water_auxiliary(soil).holding_capacity_storage,
    )
    return nothing
end

function partition_soil_water_ice!(state::ModelState)
    water_state = state.prognostic.soil.water
    water_auxiliary = state.auxiliary.soil.water
    launch_custom!(
        partition_soil_water_ice_kernel!, water_state.storage,
        size(water_state.storage, 2), water_state.ice_storage,
        water_state.wilting_ice_fraction, water_state.available_ice_storage,
        water_state.free_ice_storage, water_auxiliary.relative_content,
        water_auxiliary.free_water, water_auxiliary.wilting_storage,
        water_auxiliary.holding_capacity_storage,
    )
    return nothing
end

@kernel inbounds = true function partition_soil_water_ice_kernel!(
    liquid_storage::AbstractArray{T},
    total_ice::AbstractArray{T},
    pwp_ice_fraction::AbstractArray{T},
    available_ice::AbstractArray{T},
    free_ice::AbstractArray{T},
    relative_water::AbstractArray{T},
    free_water::AbstractArray{T},
    wilting_storage::AbstractArray{T},
    holding_storage::AbstractArray{T},
) where {T <: AbstractFloat}
    cell = @index(Global)
    # NeuralCrop currently uses LPJmL's fixed five hydrological soil layers.
    for layer in 1:5
        total = max(
            liquid_storage[layer, cell] + total_ice[layer, cell],
            zero(T),
        )
        ice = clamp(total_ice[layer, cell], zero(T), total)
        pwp_fraction, available_frozen, free_frozen,
        relative_liquid, gravitational_liquid = lpjml_water_ice_partition(
            total,
            ice,
            wilting_storage[layer, cell],
            holding_storage[layer, cell],
        )
        total_ice[layer, cell] = ice
        liquid_storage[layer, cell] = total - ice
        pwp_ice_fraction[layer, cell] = pwp_fraction
        available_ice[layer, cell] = available_frozen
        free_ice[layer, cell] = free_frozen
        relative_water[layer, cell] = relative_liquid
        free_water[layer, cell] = gravitational_liquid
    end
end
