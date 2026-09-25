_pathway_albedo!(::Val{:C3}, cft, crop, soil, pet, maize) =
    albedo!(cft, crop, soil, pet)
_pathway_albedo!(::Val{:C4}, cft, crop, soil, pet, maize) =
    albedo!(cft, crop, soil, pet; maize)

_pathway_apar!(::Val{:C3}, cft, crop, pet, snow_height, maize) =
    apar_crop!(cft, crop, pet, snow_height)
_pathway_apar!(::Val{:C4}, cft, crop, pet, snow_height, maize) =
    apar_crop!(cft, crop, pet, snow_height; maize)

"""Execute one daily crop pathway over the canonical lifecycle state."""
function _daily_crop!(
    pathway::Union{Val{:C3}, Val{:C4}}, start_day, end_day,
    processes::ProcessModules, climate, state::ModelState;
    maize = true,
    irrigation = false,
    manure = false,
    fertilizer = :auto,
    with_tillage = true,
    crop_resp_fix = false,
    nitrogen_limit_vcmax = false,
    sowing_mode::Symbol = :prescribed_sdate,
    update_vernalization_requirement::Bool = true,
    water_balance = nothing,
    nitrogen_balance = nothing,
    carbon_balance = nothing,
    thermal_balance = nothing,
    simulation_day_offset::Integer = 0,
    diagnostic_offset::Integer = 0,
    reuse_output::Bool = false,
    selected_output::Union{Nothing, Set{Tuple{Symbol, Symbol}}} = nothing,
    prescribed_phu = nothing,
    prescribed_winter_type = nothing,
    c_shift_response_sum = nothing,
    neural_parameters = nothing,
    neural_layout = nothing,
)
    cftparameters = processes.crop
    model_parameters = processes.global_parameters
    climbuf = state.prognostic.climate
    pet = state.auxiliary.pet
    managed_land = state.inputs.management
    dailyWeather = state.inputs.weather
    output = state.output

    T = eltype(crop_prognostic(state).canopy.lai)
    cftparameters = convert_precision(T, cftparameters)
    model_parameters = convert_precision(T, model_parameters)
    global_params = model_parameters.lpjml
    photo_params = model_parameters.photosynthesis
    snow_params = model_parameters.snow
    thermal_params = model_parameters.soil_thermal
    decomp_params = model_parameters.soil_decomposition
    fertilizer = fertilizer_mode(fertilizer)
    automatic_fertilizer = fertilizer === :auto
    neural_enabled = !isnothing(neural_parameters)
    neural_enabled == !isnothing(neural_layout) || throw(ArgumentError(
        "neural_parameters and neural_layout must be provided together",
    ))
    if neural_enabled
        length(neural_parameters) == neural_parameter_count(neural_layout) ||
            throw(DimensionMismatch("neural parameters do not match NeuralCrop layout"))
        components = neural_layout.components
        (!components.lambda && !components.vcmax && !components.allocation &&
         !components.decomposition && !components.snowmelt && !components.evaporation) ||
            throw(ArgumentError(
                "the production lifecycle supports only GPP, respiration, and " *
                "transpiration neural components",
            ))
    end

    if water_balance !== nothing && irrigation
        throw(ArgumentError("water-balance diagnostics currently support rainfed simulations only"))
    end

    annual_rows = count(
        climate_day -> (climate_day + simulation_day_offset) % 365 == 0,
        start_day:end_day,
    )
    output_rows = prepare_output_block!(
        output, end_day - start_day + 1, annual_rows;
        reuse = reuse_output, selected = selected_output,
    )
    annual_output_offset = 0

    for climate_day in start_day:end_day
        day = climate_day + simulation_day_offset
        block_day = climate_day - start_day + 1
        diagnostic_day = diagnostic_offset + block_day
        output_row = output_rows.first_daily_row + block_day - 1
        day_of_year = day % 365 != 0 ? day % 365 : 365
        current_co2 = readclimate!(climate, dailyWeather, climate_day)

        if carbon_balance !== nothing
            record_carbon_balance_start!(carbon_balance, diagnostic_day, state, state)
        end
        if nitrogen_balance !== nothing
            record_nitrogen_balance_start!(nitrogen_balance, diagnostic_day, state, state)
        end
        if water_balance !== nothing
            record_water_balance_start!(
                water_balance, diagnostic_day, state, dailyWeather.prec,
            )
        end

        # --- Discrete establishment event ---------------------------------
        # Today's climate must enter history before sowing decisions. This is
        # intentionally separate from continuous crop/soil process kernels.
        dynamic_winter_type = isnothing(prescribed_winter_type) ?
            crop_phenology_input(state).winter_type : prescribed_winter_type
        update_climbuf!(
            cftparameters, dailyWeather.temp, climbuf, day;
            prec = dailyWeather.prec,
            dynamic_sowing = sowing_mode === :dynamic_sdate,
            winter_type = dynamic_winter_type,
            update_vernalization_requirement,
        )
        if sowing_mode === :dynamic_sdate
            dynamic_sowing_date!(
                state, climbuf, cftparameters, day_of_year;
                irrigated = irrigation,
                prescribed_winter_type,
            )
        end
        cultivate!(
            state, managed_land, state, day_of_year;
            manure,
            apply_prescribed_fertilizer = fertilizer === :yes,
            prescribed_phu,
            prescribed_winter_type,
            cftparameters = cftparameters,
            lpjmlparams = global_params,
            laimax = cftparameters.laimax,
        )
        if carbon_balance !== nothing
            record_carbon_balance_after_cultivate!(carbon_balance, diagnostic_day, state)
        end

        if with_tillage
            litter_tillage!(state, state)
            tillage_hydraulics!(state, state; lpjmlparams = global_params)
        end
        litter_bioturbation!(state; lpjmlparams = global_params)

        # Radiation uses the snow state present at the start of the day.
        _pathway_albedo!(pathway, cftparameters, state, state, pet, maize)
        petpar!(
            pet, day_of_year, managed_land.latitude, dailyWeather.temp,
            dailyWeather.lwr, dailyWeather.swr,
        )
        sowing_mode === :dynamic_sdate &&
            record_potential_evaporation!(climbuf, pet.eeq, day, global_params)
        snow!(state, dailyWeather; snowparams = snow_params, lpjmlparams = global_params)
        if water_balance !== nothing
            record_water_balance_after_snow!(
                water_balance, diagnostic_day, dailyWeather.prec,
            )
        end

        pedotransfer!(state; lpjmlparams = global_params)
        update_surface_litter_properties!(state; thermalparams = thermal_params)
        soil_temperature!(
            state, dailyWeather.temp, climbuf.atemp_mean;
            thermalparams = thermal_params, snowparams = snow_params,
        )

        # Existing litter/SOM decomposes before today's crop uptake.
        soil_cn_decomposition!(
            state;
            lpjmlparams = global_params,
            soil_decomp_params = decomp_params,
        )
        # Match LPJmL: add atmospheric mineral N after litter/SOM turnover and
        # before the crop's daily nitrogen uptake and nitrogen losses.
        if hasproperty(climate, :no3_deposition) || hasproperty(climate, :nh4_deposition)
            nitrogen_deposition!(
                state, dailyWeather.no3_deposition, dailyWeather.nh4_deposition,
            )
        end
        c_shift_response_sum === nothing ||
            accumulate_c_shift_response!(c_shift_response_sum, state)

        # --- Discrete calendar harvest event ------------------------------
        phenology_crop!(
            state, climbuf.V_req, cftparameters, dailyWeather.temp, pet.daylength,
        )
        annual_output_row = day_of_year == 365 ?
            output_rows.first_annual_row + annual_output_offset : nothing
        harvest_crop!(
            state, state, output, managed_land.residue_fraction, day_of_year;
            output_row, annual_output_row,
        )
        route_harvest_residues!(state, state)
        annual_output_offset += day_of_year == 365
        if carbon_balance !== nothing
            record_carbon_balance_after_harvest!(
                carbon_balance, diagnostic_day, state, state,
                managed_land.residue_fraction,
            )
        end

        interception!(
            state, cftparameters, pet.eeq, dailyWeather.prec;
            lpjmlparams = global_params,
        )
        pedotransfer!(state; lpjmlparams = global_params)
        soil_infiltration!(
            state, state, dailyWeather.prec;
            snowmelt = soil_snow_fluxes(state).melt,
            air_temperature = dailyWeather.temp,
            lpjmlparams = global_params,
            thermalparams = thermal_params,
        )
        if thermal_balance !== nothing
            record_thermal_balance!(thermal_balance, diagnostic_day, state)
        end

        _pathway_apar!(
            pathway, cftparameters, state, pet,
            soil_snow_prognostic(state).height, maize,
        )
        temp_stress(
            cftparameters, pet, state, dailyWeather.temp;
            photoparams = photo_params,
        )
        photosynthesis!(
            pathway, cftparameters, state, crop_canopy_auxiliary(state).apar,
            pet.daylength, dailyWeather.temp, current_co2;
            comp_vcmax = true,
            lpjmlparams = global_params,
            photoparams = photo_params,
        )
        direct_neural_gpp = neural_enabled && neural_layout.components.gpp &&
            !neural_layout.components.gpp_residual
        direct_neural_gpp && neural_gpp!(
            neural_parameters,
            neural_layout,
            state,
            pet.daylength,
            dailyWeather.temp,
            current_co2,
            soil_properties(state).layer_depth,
            photo_params,
        )

        if neural_enabled && neural_layout.components.transpiration
            neural_transpiration!(
                neural_parameters,
                neural_layout,
                crop_fluxes(state).carbon.water_limited_assimilation,
                cftparameters,
                state,
                pet,
                current_co2,
                global_params,
            )
        else
            transpiration!(
                crop_fluxes(state).carbon.water_limited_assimilation,
                cftparameters, state, pet, state, current_co2;
                lpjmlparams = global_params,
            )
        end
        if !direct_neural_gpp
            solve_lambda!(
                pathway, cftparameters, state, pet, dailyWeather.temp, current_co2;
                lpjmlparams = global_params,
                photoparams = photo_params,
            )
        end

        if nitrogen_limit_vcmax
            crop_nitrogen!(
                state, cftparameters, state,
                crop_photosynthesis_auxiliary(state).potential_vcmax,
                dailyWeather.temp;
                auto_fertilizer = automatic_fertilizer,
                lpjmlparams = global_params,
            )
            limit_vcmax_by_nitrogen!(
                state, cftparameters, dailyWeather.temp;
                lpjmlparams = global_params,
            )
        end
        if neural_enabled && neural_layout.components.gpp
            if neural_layout.components.gpp_residual || nitrogen_limit_vcmax
                photosynthesis!(
                    pathway, cftparameters, state, crop_canopy_auxiliary(state).apar,
                    pet.daylength, dailyWeather.temp, current_co2;
                    comp_vcmax = false,
                    lpjmlparams = global_params,
                    photoparams = photo_params,
                )
                neural_gpp!(
                    neural_parameters,
                    neural_layout,
                    state,
                    pet.daylength,
                    dailyWeather.temp,
                    current_co2,
                    soil_properties(state).layer_depth,
                    photo_params,
                )
            end
        else
            photosynthesis!(
                pathway, cftparameters, state, crop_canopy_auxiliary(state).apar,
                pet.daylength, dailyWeather.temp, current_co2;
                comp_vcmax = false,
                lpjmlparams = global_params,
                photoparams = photo_params,
            )
        end

        if neural_enabled && neural_layout.components.respiration
            neural_crop_carbon!(
                neural_parameters,
                neural_layout,
                state,
                output,
                cftparameters,
                dailyWeather.temp,
                soil_thermal_prognostic(state).temperature;
                output_row,
                crop_resp_fix,
                lpjmlparams = global_params,
            )
        else
            crop_carbon!(
                state, output, cftparameters, dailyWeather.temp,
                soil_thermal_prognostic(state).temperature;
                output_row, crop_resp_fix, lpjmlparams = global_params,
            )
        end
        # --- Discrete failed-crop termination event -----------------------
        # This remains separate from calendar harvest: it is only triggered
        # after today's carbon allocation identifies a non-viable stand.
        terminate_failed_crop!(
            state, state, output, managed_land.residue_fraction, day_of_year;
            output_row, annual_output_row,
        )
        route_harvest_residues!(state, state)
        if carbon_balance !== nothing
            record_carbon_balance_after_harvest!(
                carbon_balance, diagnostic_day, state, state,
                managed_land.residue_fraction,
            )
        end

        if nitrogen_limit_vcmax
            allocate_crop_nitrogen!(state, cftparameters)
        else
            crop_nitrogen!(
                state, cftparameters, state,
                crop_photosynthesis_auxiliary(state).vcmax, dailyWeather.temp;
                auto_fertilizer = automatic_fertilizer,
                lpjmlparams = global_params,
            )
        end

        evaporation!(pet.eeq, state, state; lpjmlparams = global_params)
        soil_evapotranspiration!(state, state; irrigation)
        record_ecosystem_flux_outputs!(output, state, state; output_row)
        accumulate_season_process_diagnostics!(output, state, state)
        post_crop_nitrogen_losses!(
            state;
            air_temperature = dailyWeather.temp,
            wind_speed = dailyWeather.wind,
            lpjmlparams = global_params,
        )

        if water_balance !== nothing
            record_water_balance_end!(water_balance, diagnostic_day, state, state)
        end
        if nitrogen_balance !== nothing
            record_nitrogen_balance_end!(
                nitrogen_balance, diagnostic_day, state, state;
                nitrate_deposition = dailyWeather.no3_deposition,
                ammonium_deposition = dailyWeather.nh4_deposition,
            )
        end
        if carbon_balance !== nothing
            record_carbon_balance_end!(carbon_balance, diagnostic_day, state, state)
        end
        # All process kernels for a day share the current backend stream. Wait
        # once at the lifecycle boundary so callers observe a completed daily
        # transition without forcing a host/device barrier after every process.
        synchronize_backend!(crop_prognostic(state).canopy.lai)
    end
    return nothing
end

"""Compatibility entry point for the C3-specialized daily crop transition."""
daily_crop_C3!(args...; kwargs...) = _daily_crop!(Val(:C3), args...; kwargs...)

"""Compatibility entry point for the C4-specialized daily crop transition."""
daily_crop_C4!(args...; kwargs...) = _daily_crop!(Val(:C4), args...; kwargs...)
