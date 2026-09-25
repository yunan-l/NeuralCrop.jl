@inline function _neural_top3_water_fraction(storage, layer_depth, cell::Integer)
    T = eltype(storage)
    top3_depth = T(layer_depth[1]) + T(layer_depth[2]) + T(layer_depth[3])
    return clamp(
        (storage[1, cell] + storage[2, cell] + storage[3, cell]) /
            max(T(top3_depth), eps(T)),
        zero(T),
        one(T),
    )
end

@kernel inbounds = true function neural_set_gpp_kernel!(
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
    layout = NeuralCrop.kernel_value(fixed_layout)
    layer_depth = NeuralCrop.kernel_value(fixed_layer_depth)
    photo_parameters = NeuralCrop.kernel_value(fixed_photo_parameters)
    co2_cell = T(co2[length(co2) == 1 ? 1 : cell])
    process_gpp = T(gross_assimilation[cell])
    active = is_growing[cell] == 1 &&
        temperature_stress[cell] >= T(1e-2) &&
        daylength[cell] > zero(T) && apar[cell] > zero(T)
    gpp = if !active
        zero(T)
    elseif layout.components.gpp_residual
        process_gpp * NeuralCrop.neural_gpp_multiplier(
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
        NeuralCrop.neural_gpp(
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
    net_assimilation[cell], daily_net = NeuralCrop.compute_net_assimilation(
        gpp, T(leaf_respiration[cell]), T(daylength[cell]),
    )
    water_limited_assimilation[cell] = NeuralCrop.compute_water_limited_assimilation(
        daily_net,
        T(photo_parameters.cmass),
        T(air_temperature[cell]),
        T(photo_parameters.p),
    )
end

function _neural_set_gpp!(
    theta::AbstractVector{T},
    layout::NeuralCrop.NeuralCropLayout,
    state::NeuralCrop.ModelState,
    daylength,
    air_temperature,
    co2,
    layer_depth,
    photo_parameters,
) where {T <: AbstractFloat}
    fluxes = NeuralCrop.crop_fluxes(state).carbon
    canopy = NeuralCrop.crop_canopy_auxiliary(state)
    crop = NeuralCrop.crop_prognostic(state)
    photosynthesis = NeuralCrop.crop_photosynthesis_auxiliary(state)
    NeuralCrop.launch_1D!(
        neural_set_gpp_kernel!,
        fluxes.gross_assimilation,
        fluxes.net_assimilation,
        fluxes.water_limited_assimilation,
        theta,
        NeuralCrop.kernel_constant(theta, layout),
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
        NeuralCrop.soil_water_prognostic(state).storage,
        NeuralCrop.kernel_constant(layer_depth),
        NeuralCrop.kernel_constant(photo_parameters),
    )
    return nothing
end

@kernel inbounds = true function neural_set_lambda_kernel!(
    lambda,
    theta,
    fixed_layout,
    temperature_stress,
    apar,
    is_growing,
    daylength,
    air_temperature,
    storage,
    fixed_layer_depth,
)
    cell = @index(Global, Linear)
    T = eltype(lambda)
    layout = NeuralCrop.kernel_value(fixed_layout)
    layer_depth = NeuralCrop.kernel_value(fixed_layer_depth)
    active = is_growing[cell] == 1 &&
        temperature_stress[cell] >= T(1e-2) &&
        daylength[cell] > zero(T) && apar[cell] > zero(T)
    lambda[cell] = active ? NeuralCrop.neural_lambda(
        theta,
        layout,
        T(daylength[cell]),
        T(air_temperature[cell]),
        _neural_top3_water_fraction(storage, layer_depth, cell),
    ) : zero(T)
end

function _neural_set_lambda!(
    theta::AbstractVector{T},
    layout::NeuralCrop.NeuralCropLayout,
    state::NeuralCrop.ModelState,
    daylength,
    air_temperature,
    layer_depth,
) where {T <: AbstractFloat}
    photosynthesis = NeuralCrop.crop_photosynthesis_auxiliary(state)
    apar = NeuralCrop.crop_canopy_auxiliary(state).apar
    is_growing = NeuralCrop.crop_prognostic(state).phenology.is_growing
    storage = NeuralCrop.soil_water_prognostic(state).storage
    NeuralCrop.launch_1D!(
        neural_set_lambda_kernel!,
        photosynthesis.lambda,
        theta,
        NeuralCrop.kernel_constant(theta, layout),
        photosynthesis.temperature_stress,
        apar,
        is_growing,
        daylength,
        air_temperature,
        storage,
        NeuralCrop.kernel_constant(layer_depth),
    )
    return nothing
end

@kernel inbounds = true function neural_set_vcmax_kernel!(
    vcmax,
    theta,
    fixed_layout,
    potential_vcmax,
    nitrogen_limitation,
    lambda,
    temperature_stress,
    apar,
    leaf_nitrogen,
    is_growing,
    daylength,
    air_temperature,
    fixed_parameters,
)
    cell = @index(Global, Linear)
    T = eltype(vcmax)
    layout = NeuralCrop.kernel_value(fixed_layout)
    lpjml_params = NeuralCrop.kernel_value(fixed_parameters)
    active = is_growing[cell] == 1 &&
        temperature_stress[cell] >= T(1e-2) &&
        daylength[cell] > zero(T) && apar[cell] > zero(T)
    value = active ? NeuralCrop.neural_vcmax(
        theta,
        layout,
        T(daylength[cell]),
        T(apar[cell]),
        T(leaf_nitrogen[cell]),
        T(air_temperature[cell]),
    ) : zero(T)
    potential_vcmax[cell] = value
    vcmax[cell] = value
    nitrogen_limitation[cell] = value > zero(T) ? one(T) : zero(T)
    lambda[cell] = active ? T(lpjml_params.LAMBDA_OPT) : zero(T)
end

function _neural_set_vcmax!(
    theta::AbstractVector{T},
    layout::NeuralCrop.NeuralCropLayout,
    state::NeuralCrop.ModelState,
    daylength,
    air_temperature,
    lpjml_params,
) where {T <: AbstractFloat}
    photosynthesis = NeuralCrop.crop_photosynthesis_auxiliary(state)
    apar = NeuralCrop.crop_canopy_auxiliary(state).apar
    leaf_nitrogen = NeuralCrop.crop_prognostic(state).nitrogen.leaf
    is_growing = NeuralCrop.crop_prognostic(state).phenology.is_growing
    NeuralCrop.launch_1D!(
        neural_set_vcmax_kernel!,
        photosynthesis.vcmax,
        theta,
        NeuralCrop.kernel_constant(theta, layout),
        photosynthesis.potential_vcmax,
        photosynthesis.nitrogen_limitation,
        photosynthesis.lambda,
        photosynthesis.temperature_stress,
        apar,
        leaf_nitrogen,
        is_growing,
        daylength,
        air_temperature,
        NeuralCrop.kernel_constant(lpjml_params),
    )
    return nothing
end

@kernel inbounds = true function neural_snow_kernel!(
    pack,
    theta,
    fixed_layout,
    precipitation,
    temperature,
    melt_flux,
    sublimation_flux,
    runoff_flux,
    height,
    fraction,
    fixed_parameters,
)
    cell = @index(Global, Linear)
    T = eltype(pack)
    layout = NeuralCrop.kernel_value(fixed_layout)
    parameters = NeuralCrop.kernel_value(fixed_parameters)
    snow_params = parameters.snow_params
    lpjml_params = parameters.lpjml_params
    melt_flux[cell] = zero(T)
    sublimation_flux[cell] = zero(T)
    runoff_flux[cell] = zero(T)
    updated_pack, liquid_precipitation, runoff = NeuralCrop.compute_snowfall(
        T(pack[cell]), T(precipitation[cell]), T(temperature[cell]),
        T(snow_params.tsnow), T(lpjml_params.maxsnowpack),
    )
    updated_pack, sublimation = NeuralCrop.compute_snow_sublimation(updated_pack)
    melt = NeuralCrop.neural_snow_melt(
        theta, layout, T(temperature[cell]), updated_pack, liquid_precipitation,
    )
    updated_pack -= melt
    precipitation[cell] = liquid_precipitation + melt
    updated_pack, snow_height, snow_fraction = NeuralCrop.compute_snow_geometry(
        updated_pack, T(snow_params.c_watertosnow), T(snow_params.c_roughness),
    )
    pack[cell] = updated_pack
    height[cell] = snow_height
    fraction[cell] = snow_fraction
    melt_flux[cell] = melt
    sublimation_flux[cell] = sublimation
    runoff_flux[cell] = runoff
end

function _neural_snow!(
    theta::AbstractVector{T},
    layout::NeuralCrop.NeuralCropLayout,
    state::NeuralCrop.ModelState,
    snow_params,
    lpjml_params,
) where {T <: AbstractFloat}
    weather = state.inputs.weather
    snow = NeuralCrop.soil_snow_prognostic(state)
    flux = NeuralCrop.soil_snow_fluxes(state)
    NeuralCrop.launch_1D!(
        neural_snow_kernel!,
        snow.pack,
        theta,
        NeuralCrop.kernel_constant(theta, layout),
        weather.prec,
        weather.temp,
        flux.melt,
        flux.sublimation,
        flux.runoff,
        snow.height,
        snow.fraction,
        NeuralCrop.kernel_constant((; snow_params, lpjml_params)),
    )
    return nothing
end

@kernel inbounds = true function neural_transpiration_kernel!(
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
    layout = NeuralCrop.kernel_value(fixed_layout)
    parameters = NeuralCrop.kernel_value(fixed_parameters)
    cft = parameters.cft
    lpjml_params = parameters.lpjml_params
    soil_layers = parameters.soil_layers
    co2_index = length(co2) == 1 ? 1 : cell
    conductance = NeuralCrop.compute_canopy_conductance(
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
        supply = NeuralCrop.compute_transpiration_supply(
            T(cft.emax), root_water, root_carbon[cell],
        )
        demand = NeuralCrop.compute_transpiration_demand(
            canopy_wet[cell], equilibrium_pet[cell], T(lpjml_params.ALPHAM),
            T(lpjml_params.GM), conductance,
        )
        multiplier = NeuralCrop.neural_transpiration_multiplier(
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
        water_deficit[cell] = NeuralCrop.compute_water_sufficiency(
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
                layer_flux, capped = NeuralCrop.compute_layer_transpiration(
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
        canopy_conductance[cell] = NeuralCrop.compute_actual_canopy_conductance(
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
            transpiration_layer[layer, cell], _ = NeuralCrop.compute_layer_transpiration(
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

function _neural_transpiration!(
    theta::AbstractVector{T},
    layout::NeuralCrop.NeuralCropLayout,
    water_limited_assimilation,
    cft,
    state::NeuralCrop.ModelState,
    pet,
    co2,
    lpjml_params,
) where {T <: AbstractFloat}
    canopy = NeuralCrop.crop_canopy_auxiliary(state)
    crop = NeuralCrop.crop_prognostic(state)
    crop_water = crop.water
    crop_water_flux = NeuralCrop.crop_fluxes(state).water
    soil_water = NeuralCrop.soil_water_auxiliary(state)
    NeuralCrop.launch_1D!(
        neural_transpiration_kernel!,
        canopy.canopy_conductance,
        theta,
        NeuralCrop.kernel_constant(theta, layout),
        water_limited_assimilation,
        co2,
        pet.daylength,
        canopy.fpar,
        crop_water_flux.transpiration_layer,
        crop_water.demand_sum,
        crop_water.supply_sum,
        NeuralCrop.crop_stress_auxiliary(state).water_deficit,
        crop_water.sufficiency,
        crop.carbon.root,
        canopy.canopy_wet,
        crop.phenology.is_growing,
        pet.eeq,
        NeuralCrop.crop_root_input(state).distribution,
        NeuralCrop.crop_root_auxiliary(state).zone_available_water,
        soil_water.relative_content,
        soil_water.holding_capacity_storage,
        NeuralCrop.kernel_constant((; cft, lpjml_params, soil_layers = 5)),
    )
    return nothing
end

@kernel inbounds = true function neural_crop_respiration_kernel!(
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
    is_growing,
)
    cell = @index(Global, Linear)
    T = eltype(respiration)
    layout = NeuralCrop.kernel_value(fixed_layout)
    multiplier = NeuralCrop.neural_crop_respiration_multiplier(
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

function _neural_crop_carbon!(
    theta::AbstractVector{T},
    layout::NeuralCrop.NeuralCropLayout,
    state::NeuralCrop.ModelState,
    cft,
    air_temperature,
    soil_temperature,
    lpjml_params,
) where {T <: AbstractFloat}
    fluxes = NeuralCrop.crop_fluxes(state)
    crop = NeuralCrop.crop_prognostic(state)
    NeuralCrop.respiration!(
        state,
        cft,
        air_temperature,
        soil_temperature,
        fluxes.carbon.gross_assimilation,
        fluxes.carbon.leaf_respiration;
        lpjmlparams = lpjml_params,
    )
    NeuralCrop.launch_1D!(
        neural_crop_respiration_kernel!,
        fluxes.carbon.respiration,
        theta,
        NeuralCrop.kernel_constant(theta, layout),
        fluxes.carbon.gross_assimilation,
        fluxes.carbon.leaf_respiration,
        air_temperature,
        soil_temperature,
        crop.carbon.root,
        crop.carbon.storage,
        crop.carbon.pool,
        crop.phenology.is_growing,
    )
    NeuralCrop.carbon_allocation!(cft, state)
    return nothing
end

@kernel inbounds = true function neural_decomposition_kernel!(
    response,
    theta,
    fixed_layout,
    litter_response,
    workspace,
    water,
    thermal_surface,
    fixed_parameters,
)
    cell = @index(Global, Linear)
    T = eltype(response)
    layout = NeuralCrop.kernel_value(fixed_layout)
    decomposition_params = NeuralCrop.kernel_value(fixed_parameters)
    epsilon = T(decomposition_params.eps)
    top_moisture = zero(T)
    top_temperature = T(thermal_surface.soil_temperature[1, cell])
    for layer in axes(response, 1)
        wilting_ice = water.wilting_storage[layer, cell] *
            water.wilting_ice_fraction[layer, cell]
        liquid_capacity = water.saturation_storage[layer, cell] - wilting_ice -
            water.available_ice_storage[layer, cell] - water.free_ice_storage[layer, cell]
        moisture = (
            water.relative_content[layer, cell] * water.holding_capacity_storage[layer, cell] +
            water.wilting_storage[layer, cell] - wilting_ice + water.free_water[layer, cell]
        ) / max(liquid_capacity, T(epsilon))
        moisture = clamp(moisture, T(epsilon), one(T))
        temperature = T(thermal_surface.soil_temperature[layer, cell])
        workspace.layer_moisture[layer, cell] = moisture
        workspace.layer_temperature[layer, cell] = temperature
        response[layer, cell] = NeuralCrop.neural_decomposition_response(
            theta, layout, temperature, moisture,
        )
        if layer == 1
            top_moisture = moisture
            top_temperature = temperature
        end
    end

    surface_moisture = thermal_surface.surface_water_capacity[cell] > T(epsilon) ?
        clamp(
            thermal_surface.surface_water_storage[cell] /
                max(thermal_surface.surface_water_capacity[cell], T(epsilon)),
            zero(T),
            one(T),
        ) : top_moisture
    surface_temp = T(thermal_surface.surface_temperature[cell])
    surface_response = NeuralCrop.neural_decomposition_response(
        theta, layout, surface_temp, surface_moisture,
    )
    workspace.surface_moisture[cell] = surface_moisture
    workspace.surface_response[cell] = surface_response
    litter_response[1, cell] = surface_response
    top_response = NeuralCrop.neural_decomposition_response(
        theta, layout, top_temperature, top_moisture,
    )
    litter_response[2, cell] = top_response
    litter_response[3, cell] = top_response
end

function _neural_soil_decomposition_response!(
    theta::AbstractVector{T},
    layout::NeuralCrop.NeuralCropLayout,
    state::NeuralCrop.ModelState,
    decomposition_params,
) where {T <: AbstractFloat}
    response = NeuralCrop.soil_decomposition_auxiliary(state).response
    litter_response = NeuralCrop.soil_decomposition_auxiliary(state).litter_response
    workspace = NeuralCrop.soil_decomposition_workspace(state)
    water = NeuralCrop.soil_water_prognostic(state)
    water_auxiliary = NeuralCrop.soil_water_auxiliary(state)
    soil_temperature = NeuralCrop.soil_thermal_prognostic(state).temperature
    surface = NeuralCrop.soil_surface_litter_prognostic(state)
    surface_auxiliary = NeuralCrop.soil_surface_litter_auxiliary(state)
    NeuralCrop.launch_custom!(
        neural_decomposition_kernel!,
        response,
        size(response, 2),
        theta,
        NeuralCrop.kernel_constant(theta, layout),
        litter_response,
        (
            layer_moisture = workspace.layer_scratch_1,
            layer_temperature = workspace.layer_scratch_2,
            surface_moisture = workspace.surface_scratch_1,
            surface_response = workspace.surface_scratch_2,
        ),
        (
            relative_content = water_auxiliary.relative_content,
            saturation_storage = water_auxiliary.saturation_storage,
            holding_capacity_storage = water_auxiliary.holding_capacity_storage,
            wilting_storage = water_auxiliary.wilting_storage,
            free_water = water_auxiliary.free_water,
            wilting_ice_fraction = water.wilting_ice_fraction,
            available_ice_storage = water.available_ice_storage,
            free_ice_storage = water.free_ice_storage,
        ),
        (
            soil_temperature,
            surface_water_storage = surface.water_storage,
            surface_temperature = surface.temperature,
            surface_water_capacity = surface_auxiliary.water_capacity,
        ),
        NeuralCrop.kernel_constant(decomposition_params),
    )
    return nothing
end

@kernel inbounds = true function neural_evaporation_kernel!(
    evaporation,
    theta,
    fixed_layout,
    pet_eeq,
    fpar,
    canopy_wet,
    litter_cover,
    litter_evaporation,
    transpiration,
    relative_content,
    holding_capacity_storage,
    free_water,
    fixed_parameters,
)
    cell = @index(Global, Linear)
    T = eltype(evaporation)
    layout = NeuralCrop.kernel_value(fixed_layout)
    lpjml_params = NeuralCrop.kernel_value(fixed_parameters)
    priestley_taylor = T(lpjml_params.PRIESTLEY_TAYLOR)
    liquid_1 = relative_content[1, cell] * holding_capacity_storage[1, cell] +
        free_water[1, cell]
    liquid_2 = relative_content[2, cell] * holding_capacity_storage[2, cell] +
        free_water[2, cell]
    liquid_3 = relative_content[3, cell] * holding_capacity_storage[3, cell] +
        free_water[3, cell]
    available_1 = max(liquid_1 - transpiration[1, cell], zero(T))
    available_2 = max(liquid_2 - transpiration[2, cell], zero(T))
    available_3 = max(liquid_3 - transpiration[3, cell], zero(T))
    output = NeuralCrop.neural_soil_evaporation(
        theta,
        layout,
        T(pet_eeq[cell]),
        T(fpar[cell]),
        T(canopy_wet[cell]),
        T(litter_cover[cell]),
        T(liquid_1),
        T(liquid_2),
        T(liquid_3),
    )
    transpiration_total = zero(T)
    for layer in axes(transpiration, 1)
        transpiration_total += transpiration[layer, cell]
    end
    available_energy = max(
        T(pet_eeq[cell]) * T(priestley_taylor) *
            (one(T) - canopy_wet[cell]) - transpiration_total -
            litter_evaporation[cell],
        zero(T),
    )
    bare_soil_energy = T(pet_eeq[cell]) * T(priestley_taylor) *
        max(one(T) - fpar[cell], T(0.05))
    weighted_available = output[2] * available_1 +
        output[3] * available_2 + output[4] * available_3
    total = output[1] * min(available_energy, bare_soil_energy, weighted_available)
    if weighted_available > eps(T)
        evaporation[1, cell] = total * output[2] * available_1 / weighted_available
        evaporation[2, cell] = total * output[3] * available_2 / weighted_available
        evaporation[3, cell] = total * output[4] * available_3 / weighted_available
    else
        evaporation[1, cell] = zero(T)
        evaporation[2, cell] = zero(T)
        evaporation[3, cell] = zero(T)
    end
end

function _neural_evaporation!(
    theta::AbstractVector{T},
    layout::NeuralCrop.NeuralCropLayout,
    state::NeuralCrop.ModelState,
    pet_eeq,
    lpjml_params,
    layer_depth,
) where {T <: AbstractFloat}
    # Preserve litter evaporation and process-based layers four and five, then
    # replace only the first three soil-layer fluxes.
    NeuralCrop.evaporation!(pet_eeq, state, state; lpjmlparams = lpjml_params)
    canopy = NeuralCrop.crop_canopy_auxiliary(state)
    transpiration = NeuralCrop.crop_fluxes(state).water.transpiration_layer
    water_auxiliary = NeuralCrop.soil_water_auxiliary(state)
    evaporation = NeuralCrop.soil_water_fluxes(state).evaporation
    litter_flux = NeuralCrop.soil_surface_litter_fluxes(state)
    NeuralCrop.launch_custom!(
        neural_evaporation_kernel!,
        evaporation,
        size(evaporation, 2),
        theta,
        NeuralCrop.kernel_constant(theta, layout),
        pet_eeq,
        canopy.fpar,
        canopy.canopy_wet,
        NeuralCrop.soil_surface_litter_prognostic(state).cover,
        litter_flux.evaporation,
        transpiration,
        water_auxiliary.relative_content,
        water_auxiliary.holding_capacity_storage,
        water_auxiliary.free_water,
        NeuralCrop.kernel_constant(lpjml_params),
    )
    return nothing
end

function _neural_continuous_transition!(
    theta::AbstractVector{T},
    layout::NeuralCrop.NeuralCropLayout,
    state::NeuralCrop.ModelState,
    cft::NeuralCrop.CFTParameters,
    global_parameters,
    climate,
    day::Integer,
    layer_depth,
    irrigation::Bool,
    nitrogen_limit_vcmax::Bool,
    pathway::Union{Val{:C3}, Val{:C4}},
) where {T <: AbstractFloat}
    # Root distribution is static and is already initialized from this CFT in
    # the fixed-event state template. Keeping its kernel launch inside the
    # differentiated daily transition makes a constant scalar kernel argument
    # appear active to KernelAbstractions/Enzyme on CUDA.
    model_parameters = _enzyme_model_parameters(T, global_parameters)
    global_params = model_parameters.lpjml
    photo_params = model_parameters.photosynthesis
    snow_params = model_parameters.snow
    thermal_params = model_parameters.soil_thermal
    decomp_params = model_parameters.soil_decomposition
    climbuf = state.prognostic.climate
    pet = state.auxiliary.pet
    managed_land = state.inputs.management
    weather = state.inputs.weather

    current_co2 = NeuralCrop.readclimate!(climate, weather, day)
    NeuralCrop.update_climbuf!(
        cft, weather.temp, climbuf, day;
        update_vernalization_requirement = false,
    )
    NeuralCrop.litter_tillage!(state, state)
    NeuralCrop.tillage_hydraulics!(state; lpjmlparams = global_params)
    NeuralCrop.litter_bioturbation!(state; lpjmlparams = global_params)

    NeuralCrop._pathway_albedo!(pathway, cft, state, state, pet, true)
    NeuralCrop.petpar!(
        pet,
        day % 365 == 0 ? 365 : day % 365,
        managed_land.latitude,
        weather.temp,
        weather.lwr,
        weather.swr,
    )
    if layout.components.snowmelt
        _neural_snow!(theta, layout, state, snow_params, global_params)
    else
        NeuralCrop.snow!(
            state, weather; snowparams = snow_params, lpjmlparams = global_params,
        )
    end
    NeuralCrop.pedotransfer!(state; lpjmlparams = global_params)
    NeuralCrop.update_surface_litter_properties!(state; thermalparams = thermal_params)
    NeuralCrop.soil_temperature!(
        state,
        weather.temp,
        climbuf.atemp_mean;
        thermalparams = thermal_params,
        snowparams = snow_params,
    )
    if layout.components.decomposition
        _neural_soil_decomposition_response!(
            theta, layout, state, decomp_params,
        )
        NeuralCrop._soil_cn_decomposition_from_response!(
            state; lpjmlparams = global_params,
        )
    else
        NeuralCrop.soil_cn_decomposition!(
            state;
            lpjmlparams = global_params,
            soil_decomp_params = decomp_params,
        )
    end
    if hasproperty(climate, :no3_deposition) || hasproperty(climate, :nh4_deposition)
        NeuralCrop.nitrogen_deposition!(
            state, weather.no3_deposition, weather.nh4_deposition,
        )
    end
    NeuralCrop.phenology_crop!(
        state, climbuf.V_req, cft, weather.temp, pet.daylength,
    )
    NeuralCrop.interception!(
        state,
        cft,
        pet.eeq,
        weather.prec;
        lpjmlparams = global_params,
    )
    NeuralCrop.pedotransfer!(state; lpjmlparams = global_params)
    NeuralCrop.soil_infiltration!(
        state,
        state,
        weather.prec;
        snowmelt = NeuralCrop.soil_snow_fluxes(state).melt,
        air_temperature = weather.temp,
        lpjmlparams = global_params,
        thermalparams = thermal_params,
    )
    NeuralCrop._pathway_apar!(
        pathway,
        cft,
        state,
        pet,
        NeuralCrop.soil_snow_prognostic(state).height,
        true,
    )
    NeuralCrop.temp_stress(
        cft, pet, state, weather.temp; photoparams = photo_params,
    )

    if layout.components.gpp
        NeuralCrop.photosynthesis!(
            pathway,
            cft,
            state,
            NeuralCrop.crop_canopy_auxiliary(state).apar,
            pet.daylength,
            weather.temp,
            current_co2;
            comp_vcmax = true,
            lpjmlparams = global_params,
            photoparams = photo_params,
        )
        if !layout.components.gpp_residual
            _neural_set_gpp!(
                theta,
                layout,
                state,
                pet.daylength,
                weather.temp,
                current_co2,
                layer_depth,
                photo_params,
            )
        end
    elseif layout.components.vcmax
        _neural_set_vcmax!(
            theta,
            layout,
            state,
            pet.daylength,
            weather.temp,
            global_params,
        )
    else
        NeuralCrop.photosynthesis!(
            pathway,
            cft,
            state,
            NeuralCrop.crop_canopy_auxiliary(state).apar,
            pet.daylength,
            weather.temp,
            current_co2;
            comp_vcmax = true,
            lpjmlparams = global_params,
            photoparams = photo_params,
        )
    end
    if !layout.components.gpp && layout.components.lambda
        _neural_set_lambda!(
            theta, layout, state, pet.daylength, weather.temp, layer_depth,
        )
    end
    if !layout.components.gpp &&
       (layout.components.vcmax || layout.components.lambda)
        NeuralCrop.photosynthesis!(
            pathway,
            cft,
            state,
            NeuralCrop.crop_canopy_auxiliary(state).apar,
            pet.daylength,
            weather.temp,
            current_co2;
            comp_vcmax = false,
            lpjmlparams = global_params,
            photoparams = photo_params,
        )
    end
    if layout.components.transpiration
        _neural_transpiration!(
            theta,
            layout,
            NeuralCrop.crop_fluxes(state).carbon.water_limited_assimilation,
            cft,
            state,
            pet,
            current_co2,
            global_params,
        )
    else
        NeuralCrop.transpiration!(
            NeuralCrop.crop_fluxes(state).carbon.water_limited_assimilation,
            cft,
            state,
            pet,
            state,
            current_co2;
            lpjmlparams = global_params,
        )
    end
    if (!layout.components.gpp || layout.components.gpp_residual) &&
       !layout.components.lambda
        NeuralCrop.solve_lambda!(
            pathway,
            cft,
            state,
            pet,
            weather.temp,
            current_co2,
            lpjmlparams = global_params,
            photoparams = photo_params,
        )
    end
    if nitrogen_limit_vcmax
        NeuralCrop.crop_nitrogen!(
            state,
            cft,
            state,
            NeuralCrop.crop_photosynthesis_auxiliary(state).potential_vcmax,
            weather.temp;
            auto_fertilizer = false,
            lpjmlparams = global_params,
        )
        NeuralCrop.limit_vcmax_by_nitrogen!(
            state, cft, weather.temp; lpjmlparams = global_params,
        )
    end
    if layout.components.gpp
        if layout.components.gpp_residual || nitrogen_limit_vcmax
            NeuralCrop.photosynthesis!(
                pathway,
                cft,
                state,
                NeuralCrop.crop_canopy_auxiliary(state).apar,
                pet.daylength,
                weather.temp,
                current_co2;
                comp_vcmax = false,
                lpjmlparams = global_params,
                photoparams = photo_params,
            )
            _neural_set_gpp!(
                theta,
                layout,
                state,
                pet.daylength,
                weather.temp,
                current_co2,
                layer_depth,
                photo_params,
            )
        end
    else
        NeuralCrop.photosynthesis!(
            pathway,
            cft,
            state,
            NeuralCrop.crop_canopy_auxiliary(state).apar,
            pet.daylength,
            weather.temp,
            current_co2;
            comp_vcmax = false,
            lpjmlparams = global_params,
            photoparams = photo_params,
        )
    end
    # Carbon pools deliberately remain on NeuralCrop's native allocation path.
    if layout.components.respiration
        _neural_crop_carbon!(
            theta,
            layout,
            state,
            cft,
            weather.temp,
            NeuralCrop.soil_thermal_prognostic(state).temperature,
            global_params,
        )
    else
        _enzyme_crop_carbon!(
            state,
            cft,
            weather.temp,
            NeuralCrop.soil_thermal_prognostic(state).temperature,
            global_params,
        )
    end
    if nitrogen_limit_vcmax
        NeuralCrop.allocate_crop_nitrogen!(state, cft)
    else
        NeuralCrop.crop_nitrogen!(
            state,
            cft,
            state,
            NeuralCrop.crop_photosynthesis_auxiliary(state).vcmax,
            weather.temp;
            auto_fertilizer = false,
            lpjmlparams = global_params,
        )
    end
    if layout.components.evaporation
        _neural_evaporation!(
            theta, layout, state, pet.eeq, global_params, layer_depth,
        )
    else
        NeuralCrop.evaporation!(pet.eeq, state, state; lpjmlparams = global_params)
    end
    NeuralCrop.soil_evapotranspiration!(state, state; irrigation)
    NeuralCrop.post_crop_nitrogen_losses!(
        state;
        air_temperature = weather.temp,
        wind_speed = weather.wind,
        lpjmlparams = global_params,
    )
    return nothing
end

@inline function _neural_observables_at_cell(
    gross_assimilation,
    respiration,
    leaf_respiration,
    heterotrophic_respiration,
    interception,
    litter_evaporation,
    transpiration_layer,
    soil_evaporation,
    cell,
)
    gpp = gross_assimilation[cell]
    reco = respiration[cell] + leaf_respiration[cell] +
        heterotrophic_respiration[cell]
    et = interception[cell] + litter_evaporation[cell]
    for layer in axes(transpiration_layer, 1)
        et += transpiration_layer[layer, cell] + soil_evaporation[layer, cell]
    end
    return (gpp, reco, et)
end

@kernel inbounds = true function neural_loss_kernel!(
    gross_assimilation,
    loss_buffer,
    respiration,
    leaf_respiration,
    heterotrophic_respiration,
    interception,
    litter_evaporation,
    transpiration_layer,
    soil_evaporation,
    fixed_context,
    index,
)
    cell = @index(Global, Linear)
    # Station training has one independent model cell and one loss buffer per
    # cell. Keeping that dimension explicit avoids device scalar access and
    # leaves a direct path to batched independent trajectories.
    context = NeuralCrop.kernel_value(fixed_context)
    if context.growth_mask[index]
        T = eltype(loss_buffer)
        values = _neural_observables_at_cell(
            gross_assimilation,
            respiration,
            leaf_respiration,
            heterotrophic_respiration,
            interception,
            litter_evaporation,
            transpiration_layer,
            soil_evaporation,
            cell,
        )
        contribution = zero(T)
        if context.valid_masks.gpp[index]
            residual = (values[1] - T(context.observations.gpp[index])) /
                T(context.scales.gpp)
            contribution += T(context.weights.gpp) * residual * residual /
                T(context.counts.gpp)
        end
        if context.valid_masks.reco[index]
            residual = (values[2] - T(context.observations.reco[index])) /
                T(context.scales.reco)
            contribution += T(context.weights.reco) * residual * residual /
                T(context.counts.reco)
        end
        if context.valid_masks.et[index]
            residual = (values[3] - T(context.observations.et[index])) /
                T(context.scales.et)
            contribution += T(context.weights.et) * residual * residual /
                T(context.counts.et)
        end
        loss_buffer[cell] += contribution / T(context.weight_sum)
    end
end

function _neural_add_seasonal_loss!(
    loss_buffer,
    state::NeuralCrop.ModelState,
    _layer_depth,
    index::Int,
    context::NeuralCrop.NeuralSeasonContext,
)
    crop_flux = NeuralCrop.crop_fluxes(state)
    NeuralCrop.launch_1D!(
        neural_loss_kernel!,
        crop_flux.carbon.gross_assimilation,
        loss_buffer,
        crop_flux.carbon.respiration,
        crop_flux.carbon.leaf_respiration,
        NeuralCrop.soil_carbon_fluxes(state).heterotrophic_respiration,
        crop_flux.water.interception,
        NeuralCrop.soil_surface_litter_fluxes(state).evaporation,
        crop_flux.water.transpiration_layer,
        NeuralCrop.soil_water_fluxes(state).evaporation,
        NeuralCrop.kernel_constant(context),
        index,
    )
    return nothing
end

@kernel inbounds = true function neural_record_kernel!(
    gross_assimilation,
    gpp,
    reco,
    et,
    respiration,
    leaf_respiration,
    heterotrophic_respiration,
    interception,
    litter_evaporation,
    transpiration_layer,
    soil_evaporation,
    index,
)
    cell = @index(Global, Linear)
    values = _neural_observables_at_cell(
        gross_assimilation,
        respiration,
        leaf_respiration,
        heterotrophic_respiration,
        interception,
        litter_evaporation,
        transpiration_layer,
        soil_evaporation,
        cell,
    )
    # The current prediction API represents one station trajectory.
    if cell == 1
        gpp[index] = values[1]
        reco[index] = values[2]
        et[index] = values[3]
    end
end

function _neural_record_observables!(outputs, state, _layer_depth, index)
    crop_flux = NeuralCrop.crop_fluxes(state)
    NeuralCrop.launch_1D!(
        neural_record_kernel!,
        crop_flux.carbon.gross_assimilation,
        outputs.gpp,
        outputs.reco,
        outputs.et,
        crop_flux.carbon.respiration,
        crop_flux.carbon.leaf_respiration,
        NeuralCrop.soil_carbon_fluxes(state).heterotrophic_respiration,
        crop_flux.water.interception,
        NeuralCrop.soil_surface_litter_fluxes(state).evaporation,
        crop_flux.water.transpiration_layer,
        NeuralCrop.soil_water_fluxes(state).evaporation,
        index,
    )
    return nothing
end

function _neural_seasonal_loss_block_buffer!(
    loss_buffer::AbstractVector{T},
    theta::AbstractVector{T},
    state::NeuralCrop.ModelState,
    layout::NeuralCrop.NeuralCropLayout,
    cft::NeuralCrop.CFTParameters,
    global_parameters,
    climate,
    days::AbstractUnitRange{<:Integer},
    layer_depth,
    context::NeuralCrop.NeuralSeasonContext,
    day_range::UnitRange{Int},
    irrigation::Bool,
    nitrogen_limit_vcmax::Bool,
    pathway::Union{Val{:C3}, Val{:C4}},
) where {T <: AbstractFloat}
    for index in day_range
        _neural_continuous_transition!(
            theta,
            layout,
            state,
            cft,
            global_parameters,
            climate,
            days[index],
            layer_depth,
            irrigation,
            nitrogen_limit_vcmax,
            pathway,
        )
        _neural_add_seasonal_loss!(loss_buffer, state, layer_depth, index, context)
    end
    return nothing
end

function _neural_seasonal_loss_block(
    theta::AbstractVector{T},
    state::NeuralCrop.ModelState,
    layout::NeuralCrop.NeuralCropLayout,
    cft::NeuralCrop.CFTParameters,
    global_parameters,
    climate,
    days::AbstractUnitRange{<:Integer},
    layer_depth,
    context::NeuralCrop.NeuralSeasonContext,
    day_range::UnitRange{Int},
    irrigation::Bool,
    nitrogen_limit_vcmax::Bool,
    pathway::Union{Val{:C3}, Val{:C4}},
) where {T <: AbstractFloat}
    loss_buffer = zeros(T, length(NeuralCrop.crop_fluxes(state).carbon.gross_assimilation))
    _neural_seasonal_loss_block_buffer!(
        loss_buffer,
        theta,
        state,
        layout,
        cft,
        global_parameters,
        climate,
        days,
        layer_depth,
        context,
        day_range,
        irrigation,
        nitrogen_limit_vcmax,
        pathway,
    )
    NeuralCrop.synchronize_backend!(NeuralCrop.crop_fluxes(state).carbon.gross_assimilation)
    return loss_buffer[1]
end

"""Run one fixed-event NeuralCrop season and return its normalized loss."""
function NeuralCrop.enzyme_neural_seasonal_loss(
    theta::AbstractVector{T},
    state::NeuralCrop.ModelState,
    layout::NeuralCrop.NeuralCropLayout,
    cft::NeuralCrop.CFTParameters,
    global_parameters,
    climate,
    days::AbstractUnitRange{<:Integer},
    layer_depth,
    context::NeuralCrop.NeuralSeasonContext;
    irrigation::Bool = false,
    nitrogen_limit_vcmax::Bool = false,
) where {T <: AbstractFloat}
    length(theta) == NeuralCrop.neural_parameter_count(layout) || throw(DimensionMismatch(
        "theta does not match the NeuralCrop parameter layout",
    ))
    length(days) == length(context.growth_mask) || throw(DimensionMismatch(
        "days and context must have identical lengths",
    ))
    isempty(days) && throw(ArgumentError("days must not be empty"))
    return _neural_seasonal_loss_block(
        theta,
        state,
        layout,
        cft,
        global_parameters,
        climate,
        days,
        layer_depth,
        context,
        1:length(days),
        irrigation,
        nitrogen_limit_vcmax,
        _enzyme_pathway(cft),
    )
end

"""
    NeuralCrop.enzyme_neural_seasonal_loss_buffer!(loss_buffer, theta, state, ...)

Accumulate a fixed-event seasonal loss into a backend-resident buffer without
reading device memory on the host. The buffer must contain one element per
independent model cell and must be cleared by the caller before reuse. This is
the mutation-based objective used by CUDA reverse mode; the scalar-returning
API remains the CPU reference path.
"""
function NeuralCrop.enzyme_neural_seasonal_loss_buffer!(
    loss_buffer::AbstractVector{T},
    theta::AbstractVector{T},
    state::NeuralCrop.ModelState,
    layout::NeuralCrop.NeuralCropLayout,
    cft::NeuralCrop.CFTParameters,
    global_parameters,
    climate,
    days::AbstractUnitRange{<:Integer},
    layer_depth,
    context::NeuralCrop.NeuralSeasonContext;
    irrigation::Bool = false,
    nitrogen_limit_vcmax::Bool = false,
) where {T <: AbstractFloat}
    cells = length(NeuralCrop.crop_fluxes(state).carbon.gross_assimilation)
    length(loss_buffer) == cells || throw(DimensionMismatch(
        "loss_buffer must contain one value per model cell",
    ))
    length(theta) == NeuralCrop.neural_parameter_count(layout) || throw(DimensionMismatch(
        "theta does not match the NeuralCrop parameter layout",
    ))
    length(days) == length(context.growth_mask) || throw(DimensionMismatch(
        "days and context must have identical lengths",
    ))
    isempty(days) && throw(ArgumentError("days must not be empty"))
    _neural_seasonal_loss_block_buffer!(
        loss_buffer,
        theta,
        state,
        layout,
        cft,
        global_parameters,
        climate,
        days,
        layer_depth,
        context,
        1:length(days),
        irrigation,
        nitrogen_limit_vcmax,
        _enzyme_pathway(cft),
    )
    return nothing
end

function _neural_backend_loss_value(loss_buffer, reference)
    NeuralCrop.synchronize_backend!(reference)
    return only(Array(loss_buffer))
end

function _neural_seasonal_gradient_blockwise_device(
    theta::AbstractVector{T},
    state_factory::F,
    layout::NeuralCrop.NeuralCropLayout,
    cft::NeuralCrop.CFTParameters,
    global_parameters,
    climate,
    days::AbstractUnitRange{<:Integer},
    layer_depth,
    context::NeuralCrop.NeuralSeasonContext,
    ranges,
    pathway;
    irrigation::Bool,
    nitrogen_limit_vcmax::Bool,
) where {T <: AbstractFloat, F}
    state, _ = _state_and_shadow(state_factory)
    NeuralCrop.enzyme_prepare_daily_state!(state)
    reference = NeuralCrop.crop_fluxes(state).carbon.gross_assimilation
    cells = length(reference)
    snapshots = Vector{typeof(state)}(undef, length(ranges))
    forward_primal = zero(T)
    for (block_index, day_range) in enumerate(ranges)
        snapshots[block_index] = deepcopy(state)
        block_loss = similar(reference, T, cells)
        fill!(block_loss, zero(T))
        _neural_seasonal_loss_block_buffer!(
            block_loss,
            theta,
            state,
            layout,
            cft,
            global_parameters,
            climate,
            days,
            layer_depth,
            context,
            day_range,
            irrigation,
            nitrogen_limit_vcmax,
            pathway,
        )
        forward_primal += _neural_backend_loss_value(block_loss, reference)
    end

    gradient = similar(theta)
    fill!(gradient, zero(T))
    state_cotangent = nothing
    reverse_primal = zero(T)
    for block_index in length(ranges):-1:1
        block_state = deepcopy(snapshots[block_index])
        block_shadow = state_cotangent === nothing ?
            NeuralCrop.enzyme_zero_tangent(block_state) : deepcopy(state_cotangent)
        block_gradient = similar(theta)
        fill!(block_gradient, zero(T))
        block_loss = similar(reference, T, cells)
        fill!(block_loss, zero(T))
        block_loss_shadow = similar(block_loss)
        fill!(block_loss_shadow, one(T))
        # Static activity is paired with explicit KernelConstant wrappers at
        # launches: KernelAbstractions' GPU reverse rule cannot represent an
        # Active immutable scalar or parameter bundle.
        Enzyme.autodiff(
            Enzyme.Reverse,
            _neural_seasonal_loss_block_buffer!,
            Enzyme.Duplicated(block_loss, block_loss_shadow),
            Enzyme.Duplicated(theta, block_gradient),
            Enzyme.Duplicated(block_state, block_shadow),
            Enzyme.Const(layout),
            Enzyme.Const(cft),
            Enzyme.Const(global_parameters),
            Enzyme.Const(climate),
            Enzyme.Const(days),
            Enzyme.Const(layer_depth),
            Enzyme.Const(context),
            Enzyme.Const(ranges[block_index]),
            Enzyme.Const(irrigation),
            Enzyme.Const(nitrogen_limit_vcmax),
            Enzyme.Const(pathway),
        )
        reverse_primal += _neural_backend_loss_value(block_loss, reference)
        gradient .+= block_gradient
        state_cotangent = block_shadow
    end
    return (;
        primal = reverse_primal,
        forward_primal,
        gradient,
        block_days = length(first(ranges)),
        block_ranges = ranges,
    )
end

"""Compute a checkpointed Enzyme reverse gradient for all retained MLP blocks."""
function NeuralCrop.enzyme_neural_seasonal_gradient_blockwise(
    theta::AbstractVector{T},
    state_factory::F,
    layout::NeuralCrop.NeuralCropLayout,
    cft::NeuralCrop.CFTParameters,
    global_parameters,
    climate,
    days::AbstractUnitRange{<:Integer},
    layer_depth,
    context::NeuralCrop.NeuralSeasonContext;
    block_days::Integer = 30,
    irrigation::Bool = false,
    nitrogen_limit_vcmax::Bool = false,
) where {T <: AbstractFloat, F}
    length(theta) == NeuralCrop.neural_parameter_count(layout) || throw(DimensionMismatch(
        "theta does not match the NeuralCrop parameter layout",
    ))
    length(days) == length(context.growth_mask) || throw(DimensionMismatch(
        "days and context must have identical lengths",
    ))
    ranges = _seasonal_block_ranges(days, Int(block_days))
    isempty(ranges) && throw(ArgumentError("days must not be empty"))
    pathway = _enzyme_pathway(cft)

    backend = KernelAbstractions.get_backend(theta)
    if !(backend isa KernelAbstractions.CPU)
        return _neural_seasonal_gradient_blockwise_device(
            theta,
            state_factory,
            layout,
            cft,
            global_parameters,
            climate,
            days,
            layer_depth,
            context,
            ranges,
            pathway;
            irrigation,
            nitrogen_limit_vcmax,
        )
    end

    state, _ = _state_and_shadow(state_factory)
    NeuralCrop.enzyme_prepare_daily_state!(state)
    snapshots = Vector{typeof(state)}(undef, length(ranges))
    forward_primal = zero(T)
    for (block_index, day_range) in enumerate(ranges)
        snapshots[block_index] = deepcopy(state)
        forward_primal += _neural_seasonal_loss_block(
            theta,
            state,
            layout,
            cft,
            global_parameters,
            climate,
            days,
            layer_depth,
            context,
            day_range,
            irrigation,
            nitrogen_limit_vcmax,
            pathway,
        )
    end

    gradient = similar(theta)
    fill!(gradient, zero(T))
    state_cotangent = nothing
    reverse_primal = zero(T)
    for block_index in length(ranges):-1:1
        block_state = deepcopy(snapshots[block_index])
        block_shadow = state_cotangent === nothing ?
            NeuralCrop.enzyme_zero_tangent(block_state) : deepcopy(state_cotangent)
        block_gradient = similar(theta)
        fill!(block_gradient, zero(T))
        result = Enzyme.autodiff(
            Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal),
            _neural_seasonal_loss_block,
            Enzyme.Duplicated(theta, block_gradient),
            Enzyme.Duplicated(block_state, block_shadow),
            Enzyme.Const(layout),
            Enzyme.Const(cft),
            Enzyme.Const(global_parameters),
            Enzyme.Const(climate),
            Enzyme.Const(days),
            Enzyme.Const(layer_depth),
            Enzyme.Const(context),
            Enzyme.Const(ranges[block_index]),
            Enzyme.Const(irrigation),
            Enzyme.Const(nitrogen_limit_vcmax),
            Enzyme.Const(pathway),
        )
        reverse_primal += result[2]
        gradient .+= block_gradient
        state_cotangent = block_shadow
    end
    return (;
        primal = reverse_primal,
        forward_primal,
        gradient,
        block_days = Int(block_days),
        block_ranges = ranges,
    )
end

"""Run one season and return daily GPP, RECO, and ET."""
function NeuralCrop.enzyme_neural_predictions(
    theta::AbstractVector{T},
    state_factory::F,
    layout::NeuralCrop.NeuralCropLayout,
    cft::NeuralCrop.CFTParameters,
    global_parameters,
    climate,
    days::AbstractUnitRange{<:Integer},
    layer_depth;
    irrigation::Bool = false,
    nitrogen_limit_vcmax::Bool = false,
) where {T <: AbstractFloat, F}
    length(theta) == NeuralCrop.neural_parameter_count(layout) || throw(DimensionMismatch(
        "theta does not match the NeuralCrop parameter layout",
    ))
    state, _ = _state_and_shadow(state_factory)
    NeuralCrop.enzyme_prepare_daily_state!(state)
    gpp = similar(theta, T, length(days))
    reco = similar(gpp)
    et = similar(gpp)
    pathway = _enzyme_pathway(cft)
    for index in eachindex(days)
        _neural_continuous_transition!(
            theta,
            layout,
            state,
            cft,
            global_parameters,
            climate,
            days[index],
            layer_depth,
            irrigation,
            nitrogen_limit_vcmax,
            pathway,
        )
        _neural_record_observables!(
            (; gpp, reco, et), state, layer_depth, index,
        )
    end
    NeuralCrop.synchronize_backend!(NeuralCrop.crop_fluxes(state).carbon.gross_assimilation)
    return (; gpp, reco, et)
end
