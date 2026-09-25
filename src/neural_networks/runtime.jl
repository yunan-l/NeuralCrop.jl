"""Return the mean volumetric water fraction in the first three soil layers."""
@inline function _neural_top3_water_fraction(storage, layer_depth, cell::Integer)
    T = eltype(storage)
    top3_depth = T(layer_depth[1]) + T(layer_depth[2]) + T(layer_depth[3])
    return clamp(
        (storage[1, cell] + storage[2, cell] + storage[3, cell]) /
            max(top3_depth, eps(T)),
        zero(T),
        one(T),
    )
end

@kernel inbounds = true function _neural_gpp_kernel!(
    gross_assimilation,
    net_assimilation,
    water_limited_assimilation,
    theta,
    fixed_layout,
    leaf_respiration,
    temperature_stress,
    apar,
    fpar,
    lai,
    leaf_nitrogen,
    is_growing,
    daylength,
    air_temperature,
    co2,
    storage,
    fixed_layer_depth,
    fixed_photo_parameters,
)
    cell = @index(Global, Linear)
    T = eltype(gross_assimilation)
    layout = kernel_value(fixed_layout)
    layer_depth = kernel_value(fixed_layer_depth)
    photo_parameters = kernel_value(fixed_photo_parameters)
    co2_cell = T(co2[length(co2) == 1 ? 1 : cell])
    process_gpp = T(gross_assimilation[cell])
    active = is_growing[cell] == 1 &&
        temperature_stress[cell] >= T(1e-2) &&
        daylength[cell] > zero(T) && apar[cell] > zero(T)
    gpp = if !active
        zero(T)
    elseif layout.components.gpp_residual
        process_gpp * neural_gpp_multiplier(
            theta,
            layout,
            T(daylength[cell]),
            T(apar[cell]),
            T(fpar[cell]),
            T(lai[cell]),
            T(leaf_nitrogen[cell]),
            T(air_temperature[cell]),
            co2_cell,
            _neural_top3_water_fraction(storage, layer_depth, cell),
        )
    else
        neural_gpp(
            theta,
            layout,
            T(daylength[cell]),
            T(apar[cell]),
            T(fpar[cell]),
            T(lai[cell]),
            T(leaf_nitrogen[cell]),
            T(air_temperature[cell]),
            co2_cell,
            _neural_top3_water_fraction(storage, layer_depth, cell),
        )
    end
    gross_assimilation[cell] = gpp
    net_assimilation[cell], daily_net = compute_net_assimilation(
        gpp, T(leaf_respiration[cell]), T(daylength[cell]),
    )
    water_limited_assimilation[cell] = compute_water_limited_assimilation(
        daily_net,
        T(photo_parameters.cmass),
        T(air_temperature[cell]),
        T(photo_parameters.p),
    )
end

"""Replace today's process GPP while preserving downstream carbon accounting."""
function neural_gpp!(
    theta::AbstractVector{T},
    layout::NeuralCropLayout,
    state::ModelState,
    daylength,
    air_temperature,
    co2,
    layer_depth,
    photo_parameters,
) where {T <: AbstractFloat}
    fluxes = crop_fluxes(state).carbon
    canopy = crop_canopy_auxiliary(state)
    crop = crop_prognostic(state)
    photosynthesis = crop_photosynthesis_auxiliary(state)
    launch_1D!(
        _neural_gpp_kernel!,
        fluxes.gross_assimilation,
        fluxes.net_assimilation,
        fluxes.water_limited_assimilation,
        theta,
        kernel_constant(theta, layout),
        fluxes.leaf_respiration,
        photosynthesis.temperature_stress,
        canopy.apar,
        canopy.fpar,
        crop.canopy.lai,
        crop.nitrogen.leaf,
        crop.phenology.is_growing,
        daylength,
        air_temperature,
        co2,
        soil_water_prognostic(state).storage,
        kernel_constant(theta, layer_depth),
        kernel_constant(theta, photo_parameters),
    )
    return nothing
end

@kernel inbounds = true function _neural_transpiration_kernel!(
    canopy_conductance,
    theta,
    fixed_layout,
    water_limited_assimilation,
    co2,
    daylength,
    fpar,
    transpiration_layer,
    demand_sum,
    supply_sum,
    water_deficit,
    water_sufficiency,
    root_carbon,
    canopy_wet,
    is_growing,
    equilibrium_pet,
    root_distribution,
    rootzone_available_water,
    relative_water,
    holding_capacity_storage,
    fixed_parameters,
)
    cell = @index(Global, Linear)
    T = eltype(canopy_conductance)
    layout = kernel_value(fixed_layout)
    parameters = kernel_value(fixed_parameters)
    cft = parameters.cft
    lpjml_params = parameters.lpjml_params
    soil_layers = parameters.soil_layers
    co2_index = length(co2) == 1 ? 1 : cell
    conductance = compute_canopy_conductance(
        water_limited_assimilation[cell],
        co2[co2_index],
        daylength[cell],
        fpar[cell],
        T(cft.gmin),
        T(lpjml_params.LAMBDA_OPT),
    )
    canopy_conductance[cell] = conductance

    root_water = zero(T)
    observed_rootzone_water = zero(T)
    for layer in 1:soil_layers
        root_water += relative_water[layer, cell] * root_distribution[layer]
        if layer <= 3
            observed_rootzone_water += relative_water[layer, cell] *
                holding_capacity_storage[layer, cell] * root_distribution[layer]
        end
    end
    rootzone_available_water[cell] = observed_rootzone_water

    if is_growing[cell] == 1
        supply = compute_transpiration_supply(T(cft.emax), root_water, root_carbon[cell])
        demand = compute_transpiration_demand(
            canopy_wet[cell], equilibrium_pet[cell], T(lpjml_params.ALPHAM),
            T(lpjml_params.GM), conductance,
        )
        multiplier = neural_transpiration_multiplier(
            theta,
            layout,
            T(daylength[cell]),
            T(equilibrium_pet[cell]),
            T(fpar[cell]),
            T(canopy_wet[cell]),
            conductance,
            T(root_carbon[cell]),
            root_water,
        )
        potential = min(multiplier * min(supply, demand), demand)
        demand_sum[cell] += demand
        supply_sum[cell] += potential
        water_deficit[cell] = compute_water_sufficiency(
            supply_sum[cell], demand_sum[cell],
        )
        if equilibrium_pet[cell] > zero(T) && conductance > zero(T)
            base_sufficiency = T(cft.emax) * root_water /
                (equilibrium_pet[cell] * T(lpjml_params.ALPHAM) /
                 (one(T) + (T(lpjml_params.GM) * T(lpjml_params.ALPHAM)) /
                  conductance))
            water_sufficiency[cell] = clamp(
                multiplier * base_sufficiency, zero(T), one(T),
            )
        else
            water_sufficiency[cell] = one(T)
        end
        transpiration = root_water > zero(T) ?
            potential / root_water * T(cft.fpc) : zero(T)
        corrected_total = zero(T)
        if transpiration > zero(T)
            for layer in 1:soil_layers
                layer_flux, capped = compute_layer_transpiration(
                    transpiration,
                    T(root_distribution[layer]),
                    T(relative_water[layer, cell]),
                    T(holding_capacity_storage[layer, cell]),
                )
                corrected_total += layer_flux
                capped && corrected_total < T(1e-5) && (corrected_total = zero(T))
            end
        end

        actual_supply = cft.fpc > zero(T) ? corrected_total / T(cft.fpc) : zero(T)
        canopy_conductance[cell] = compute_actual_canopy_conductance(
            conductance,
            actual_supply,
            demand,
            canopy_wet[cell],
            equilibrium_pet[cell],
            T(lpjml_params.ALPHAM),
            T(lpjml_params.GM),
        )

        transpiration = root_water > zero(T) ? corrected_total / root_water : zero(T)
        for layer in 1:soil_layers
            transpiration_layer[layer, cell], _ = compute_layer_transpiration(
                transpiration,
                T(root_distribution[layer]),
                T(relative_water[layer, cell]),
                T(holding_capacity_storage[layer, cell]),
            )
        end
    else
        canopy_conductance[cell] = zero(T)
        for layer in 1:soil_layers
            transpiration_layer[layer, cell] = zero(T)
        end
        demand_sum[cell] = zero(T)
        supply_sum[cell] = zero(T)
        water_deficit[cell] = zero(T)
        water_sufficiency[cell] = one(T)
    end
end

"""Replace crop transpiration while retaining demand and layer water caps."""
function neural_transpiration!(
    theta::AbstractVector{T},
    layout::NeuralCropLayout,
    water_limited_assimilation,
    cft,
    state::ModelState,
    pet,
    co2,
    lpjml_params,
) where {T <: AbstractFloat}
    canopy = crop_canopy_auxiliary(state)
    crop = crop_prognostic(state)
    crop_water = crop.water
    crop_water_flux = crop_fluxes(state).water
    soil_water = soil_water_auxiliary(state)
    soil_layers = size(soil_water.relative_content, 1)
    launch_1D!(
        _neural_transpiration_kernel!,
        canopy.canopy_conductance,
        theta,
        kernel_constant(theta, layout),
        water_limited_assimilation,
        co2,
        pet.daylength,
        canopy.fpar,
        crop_water_flux.transpiration_layer,
        crop_water.demand_sum,
        crop_water.supply_sum,
        crop_stress_auxiliary(state).water_deficit,
        crop_water.sufficiency,
        crop.carbon.root,
        canopy.canopy_wet,
        crop.phenology.is_growing,
        pet.eeq,
        crop_root_input(state).distribution,
        crop_root_auxiliary(state).zone_available_water,
        soil_water.relative_content,
        soil_water.holding_capacity_storage,
        kernel_constant(theta, (; cft, lpjml_params, soil_layers)),
    )
    return nothing
end

@kernel inbounds = true function _neural_crop_respiration_kernel!(
    respiration,
    theta,
    fixed_layout,
    gross_assimilation,
    leaf_respiration,
    air_temperature,
    soil_temperature,
    root_carbon,
    storage_carbon,
    pool_carbon,
)
    cell = @index(Global, Linear)
    T = eltype(respiration)
    layout = kernel_value(fixed_layout)
    multiplier = neural_crop_respiration_multiplier(
        theta,
        layout,
        T(gross_assimilation[cell]),
        T(leaf_respiration[cell]),
        T(air_temperature[cell]),
        T(soil_temperature[1, cell]),
        T(root_carbon[cell]),
        T(storage_carbon[cell]),
        T(pool_carbon[cell]),
    )
    respiration[cell] *= multiplier
end

"""Replace crop respiration and keep native carbon allocation and outputs."""
function neural_crop_carbon!(
    theta::AbstractVector{T},
    layout::NeuralCropLayout,
    state::ModelState,
    output::Output,
    cft::CFTParameters,
    air_temperature::AbstractVector{T},
    soil_temperature::AbstractMatrix{T};
    output_row::Union{Nothing, Integer} = nothing,
    crop_resp_fix::Bool = false,
    lpjmlparams::LPJmLParams = lpjmlparams,
) where {T <: AbstractFloat}
    fluxes = crop_fluxes(state)
    crop = crop_prognostic(state)
    respiration!(
        state,
        cft,
        air_temperature,
        soil_temperature,
        fluxes.carbon.gross_assimilation,
        fluxes.carbon.leaf_respiration;
        crop_resp_fix,
        lpjmlparams,
    )
    launch_1D!(
        _neural_crop_respiration_kernel!,
        fluxes.carbon.respiration,
        theta,
        kernel_constant(theta, layout),
        fluxes.carbon.gross_assimilation,
        fluxes.carbon.leaf_respiration,
        air_temperature,
        soil_temperature,
        crop.carbon.root,
        crop.carbon.storage,
        crop.carbon.pool,
    )
    carbon_allocation!(cft, state)

    sources = (
        gpp = fluxes.carbon.gross_assimilation,
        npp = fluxes.carbon.npp,
        lambda = crop_photosynthesis_auxiliary(state).lambda,
        potential_vcmax = crop_photosynthesis_auxiliary(state).potential_vcmax,
        vcmax = crop_photosynthesis_auxiliary(state).vcmax,
        nitrogen_limitation = crop_photosynthesis_auxiliary(state).nitrogen_limitation,
        respiration = fluxes.carbon.respiration,
        lai = crop_canopy_auxiliary(state).actual_lai,
        fphu = crop_phenology_auxiliary(state).fphu,
        water_deficit = crop_stress_auxiliary(state).water_deficit,
        biomass = crop.carbon.biomass,
    )
    for (field, source) in pairs(sources)
        if output_row === nothing
            setproperty!(
                output.crop,
                field,
                _append_output_row(getproperty(output.crop, field), source),
            )
        else
            _write_output_row!(getproperty(output.crop, field), output_row, source)
        end
    end
    return nothing
end
