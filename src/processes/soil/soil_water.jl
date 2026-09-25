"""
    soil_infiltration!(soil, crop, precipitation)

Apply throughfall infiltration and percolation before the daily plant water-stress calculation.
Rainfed and full-irrigation simulations share this hydrologic and thermal
update. Full irrigation differs only later in `soil_evapotranspiration!`,
where storage is restored to field capacity after evapotranspiration.
"""
function soil_infiltration!(soil,
                            crop,
                            prec::AbstractArray{T};
                            snowmelt::Union{Nothing, AbstractArray{T}} = nothing,
                            air_temperature::Union{Nothing, AbstractArray{T}} = nothing,
                            lpjmlparams::LPJmLParams = lpjmlparams,
                            thermalparams::SoilThermalParams{T} = SoilThermalParams{T}(),
) where {T <: AbstractFloat}
    surface_litter_interception!(soil, prec, crop_fluxes(crop).water.interception)
    # Full irrigation is an unlimited water-supply boundary condition, not a
    # separate infiltration physics. Keep rainfall, snowmelt, percolation, and
    # their enthalpy transport identical to the rainfed path before the later
    # field-capacity restoration.
    transfer_heat = snowmelt !== nothing && air_temperature !== nothing
    if transfer_heat
        infil_perc!(
            soil, prec, snowmelt, air_temperature;
            lpjmlparams = lpjmlparams,
            thermalparams = thermalparams,
        )
    else
        infil_perc!(soil; lpjmlparams = lpjmlparams)
    end

    launch_custom!(
        add_layer_flux_kernel!,
        soil_water_prognostic(soil).storage,
        size(soil_water_prognostic(soil).storage, 2),
        soil_water_fluxes(soil).percolation,
        size(soil_water_prognostic(soil).storage, 1),
    )
    if transfer_heat
        # LPJmL may reconcile temperature every two infiltration iterations.
        # NeuralCrop preserves the same water/energy ledger but applies it once
        # after the GPU column kernel, avoiding device synchronization inside
        # the iterative hydrology loop.
        apply_percolation_enthalpy!(soil; thermalparams = thermalparams)
    else
        partition_soil_water_ice!(soil)
    end

    return nothing
end


"""
    soil_evapotranspiration!(soil, crop; irrigation=false)

Remove the current day's layer-resolved transpiration and soil evaporation after
plant water demand and supply have been calculated.
"""
function soil_evapotranspiration!(soil,
                                  crop;
                                  irrigation = false)
    if irrigation
        launch_custom!(
            reset_irrigated_storage_kernel!,
            soil_water_prognostic(soil).storage,
            size(soil_water_prognostic(soil).storage, 2),
            soil_water_prognostic(soil).ice_storage,
            soil_water_auxiliary(soil).field_capacity,
            soil_properties(soil).layer_depth,
            size(soil_water_prognostic(soil).storage, 1),
        )
    else
        launch_custom!(
            remove_evapotranspiration_kernel!,
            soil_water_prognostic(soil).storage,
            size(soil_water_prognostic(soil).storage, 2),
            crop_fluxes(crop).water.transpiration_layer,
            soil_water_fluxes(soil).evaporation,
            size(soil_water_prognostic(soil).storage, 1),
        )
    end
    partition_soil_water_ice!(soil)

    return nothing
end

@kernel inbounds = true function add_layer_flux_kernel!(
    storage::AbstractMatrix{T},
    flux::AbstractMatrix{T},
    layers::Integer,
) where {T <: AbstractFloat}
    cell = @index(Global)
    for layer in 1:layers
        storage[layer, cell] += flux[layer, cell]
    end
end

@kernel inbounds = true function remove_evapotranspiration_kernel!(
    storage::AbstractMatrix{T},
    transpiration::AbstractMatrix{T},
    evaporation::AbstractMatrix{T},
    layers::Integer,
) where {T <: AbstractFloat}
    cell = @index(Global)
    for layer in 1:layers
        storage[layer, cell] -= transpiration[layer, cell] + evaporation[layer, cell]
    end
end

@kernel inbounds = true function reset_irrigated_storage_kernel!(
    storage::AbstractMatrix{T},
    ice_storage::AbstractMatrix{T},
    field_capacity::AbstractMatrix{T},
    layer_depth::AbstractVector{T},
    layers::Integer,
) where {T <: AbstractFloat}
    cell = @index(Global)
    for layer in 1:layers
        target_storage = field_capacity[layer, cell] * layer_depth[layer]
        storage[layer, cell] = max(
            target_storage - ice_storage[layer, cell],
            zero(T),
        )
    end
end
