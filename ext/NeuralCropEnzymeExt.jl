module NeuralCropEnzymeExt

import NeuralCrop
import Enzyme
import KernelAbstractions
using KernelAbstractions: @index, @kernel

# Fixed floating-point parameter bundles otherwise become `Active` arguments
# at nested GPU kernel calls, which KernelAbstractions' Enzyme reverse rule
# cannot represent. The wrapper is opt-in at launches, so ordinary parameter
# differentiation remains available outside NeuralCrop's fixed-parameter path.
@inline Enzyme.EnzymeRules.inactive_type(
    ::Type{<:NeuralCrop.KernelConstant},
) = true

# Output buffers are diagnostics for this scalar transition objective. They
# are still executed by the primal function, but do not participate in its
# derivative. Keeping these rules in the optional extension avoids changing
# the production output path or its ordinary runtime behavior.
function Enzyme.EnzymeRules.inactive(
    ::typeof(NeuralCrop.prepare_output_block!),
    args...,
)
    return nothing
end

function Enzyme.EnzymeRules.inactive(
    ::typeof(NeuralCrop._write_output_row!),
    args...,
)
    return nothing
end

function Enzyme.EnzymeRules.inactive(
    ::typeof(NeuralCrop.record_ecosystem_flux_outputs!),
    args...,
)
    return nothing
end

function Enzyme.EnzymeRules.inactive(
    ::typeof(NeuralCrop.convert_precision),
    ::Type{<:AbstractFloat},
    ::NeuralCrop.ModelParameters,
)
    return nothing
end

# Climate forcing is a fixed input to every current Enzyme objective. Execute
# the production reader for its primal weather-buffer updates, while keeping
# its backend launch outside Enzyme's differentiated graph.
function Enzyme.EnzymeRules.inactive(::typeof(NeuralCrop.readclimate!), args...)
    return nothing
end

# Deposition is a fixed additive input to the mineral-N state for the current
# objectives. The primal kernel runs normally; its parameter tangent is the
# unchanged incoming mineral-N tangent.
function Enzyme.EnzymeRules.inactive(::typeof(NeuralCrop.nitrogen_deposition!), args...)
    return nothing
end

function _zero_tangent_arrays!(value)
    value isa AbstractArray && return fill!(value, zero(eltype(value)))
    value === nothing && return nothing
    value isa Union{Symbol, String, Function, Type} && return value
    isstructtype(typeof(value)) || return value
    for field in fieldnames(typeof(value))
        _zero_tangent_arrays!(getfield(value, field))
    end
    return value
end

"""Create an Enzyme zero tangent for a mutable NeuralCrop state container."""
function NeuralCrop.enzyme_zero_tangent(value)
    # Some nested immutable containers receive shallow array fields from
    # `make_zero`; isolate them before clearing the tangent arrays.
    tangent = deepcopy(Enzyme.make_zero(value))
    _zero_tangent_arrays!(tangent)
    return tangent
end

function _state_and_shadow(state_factory)
    pair = state_factory()
    pair isa Tuple && length(pair) == 2 || throw(ArgumentError(
        "state_factory must return (state, state_shadow)",
    ))
    state, state_shadow = pair
    typeof(state) === typeof(state_shadow) || throw(ArgumentError(
        "state and state_shadow must have identical types",
    ))
    return state, state_shadow
end

"""
    NeuralCrop.enzyme_forward_directional(
        objective, theta, direction, state_factory, args...;
        runtime_activity=true,
    )

Differentiate a scalar, in-place NeuralCrop objective in one parameter direction.
`objective(theta, state, args...)` must return a scalar and may mutate `state`.
`state_factory()` is called for every probe and must return a fresh
`(state, state_shadow)` pair. This keeps the AD call from reusing a mutated
state between parameter directions.
"""
function NeuralCrop.enzyme_forward_directional(
    objective::F,
    theta::AbstractVector{T},
    direction::AbstractVector{T},
    state_factory::SF,
    args...;
    runtime_activity::Bool = true,
    return_primal::Bool = false,
) where {F, T <: AbstractFloat, SF}
    axes(theta) == axes(direction) || throw(DimensionMismatch(
        "theta and direction must have identical axes",
    ))
    state, state_shadow = _state_and_shadow(state_factory)
    mode = if return_primal
        Enzyme.set_runtime_activity(Enzyme.ForwardWithPrimal)
    elseif runtime_activity
        Enzyme.set_runtime_activity(Enzyme.Forward)
    else
        Enzyme.Forward
    end
    annotations = map(Enzyme.Const, args)
    result = Enzyme.autodiff(
        mode,
        objective,
        Enzyme.Duplicated,
        Enzyme.Duplicated(theta, direction),
        Enzyme.Duplicated(state, state_shadow),
        annotations...,
    )
    return return_primal ? (primal = result[2], directional = result[1]) : result[1]
end

"""
    NeuralCrop.enzyme_forward_gradient(objective, theta, state_factory, args...)

Compute one forward-mode Enzyme derivative for each element of `theta`. The
state factory is called once per parameter, so each derivative starts from the
same state. This is intentionally a small validation API; reverse-mode
optimization will be added only after the in-place state contract is verified.
"""
function NeuralCrop.enzyme_forward_gradient(
    objective::F,
    theta::AbstractVector{T},
    state_factory::SF,
    args...;
    runtime_activity::Bool = true,
) where {F, T <: AbstractFloat, SF}
    gradient = similar(theta)
    direction = similar(theta)
    for index in eachindex(theta)
        fill!(direction, zero(T))
        direction[index] = one(T)
        gradient[index] = NeuralCrop.enzyme_forward_directional(
            objective,
            theta,
            direction,
            state_factory,
            args...;
            runtime_activity,
        )
    end
    return gradient
end

"""Preallocate the one-row output buffer required by the AD daily wrapper."""
function NeuralCrop.enzyme_prepare_daily_state!(state::NeuralCrop.ModelState)
    NeuralCrop.prepare_output_block!(state.output, 1, 1; reuse = true)
    return state
end

# `crop_carbon!` combines continuous carbon processes with dynamic diagnostic
# field writes. Reuse its two process calls while leaving the diagnostic path
# outside the Enzyme graph.
function _enzyme_crop_carbon!(
    state::NeuralCrop.ModelState,
    cft::NeuralCrop.CFTParameters,
    air_temperature,
    soil_temperature,
    lpjmlparams,
)
    fluxes = NeuralCrop.crop_fluxes(state)
    NeuralCrop.respiration!(
        state,
        cft,
        air_temperature,
        soil_temperature,
        fluxes.carbon.gross_assimilation,
        fluxes.carbon.leaf_respiration;
        lpjmlparams,
    )
    NeuralCrop.carbon_allocation!(cft, state)
    return nothing
end

@inline function _enzyme_c3_lambda_residual_slope(
    fac,
    lambda,
    vcmax,
    stress,
    b,
    co2,
    temperature,
    apar,
    daylength,
    lpjmlparams,
    photoparams,
)
    T = typeof(fac)
    stress < T(1e-2) && return -fac
    ko = lpjmlparams.ko25 * photoparams.q10ko ^ ((temperature - T(25)) * T(0.1))
    kc = lpjmlparams.kc25 * photoparams.q10kc ^ ((temperature - T(25)) * T(0.1))
    co2_fac = kc * (one(T) + photoparams.po2 / ko)
    tau = photoparams.tau25 * photoparams.q10tau ^ ((temperature - T(25)) * T(0.1))
    gammastar = photoparams.po2 / (T(2) * tau)
    internal_co2 = lambda * co2
    denominator_1 = internal_co2 + T(2) * gammastar
    denominator_2 = internal_co2 + co2_fac
    c1_slope = stress * lpjmlparams.alphac3 *
        co2 * T(3) * gammastar / (denominator_1 * denominator_1)
    c2_slope = co2 * (co2_fac + gammastar) / (denominator_2 * denominator_2)
    je_slope = c1_slope * apar * photoparams.cmass * photoparams.cq / daylength
    jc_slope = c2_slope * vcmax / T(24)
    je = stress * lpjmlparams.alphac3 *
        ((internal_co2 - gammastar) / denominator_1) *
        apar * photoparams.cmass * photoparams.cq / daylength
    jc = ((internal_co2 - gammastar) / denominator_2) * vcmax / T(24)
    total = je + jc
    radicand = max(zero(T), total * total - T(4) * lpjmlparams.theta * je * jc)
    root = sqrt(radicand)
    root_slope = root > zero(T) ?
        (total * (je_slope + jc_slope) -
         T(2) * lpjmlparams.theta * (je_slope * jc + je * jc_slope)) / root :
        zero(T)
    agd_slope = (je_slope + jc_slope - root_slope) /
        (T(2) * lpjmlparams.theta) * daylength
    adt_slope = agd_slope
    scale = (temperature + T(273.15)) / photoparams.p * T(8.314) / photoparams.cmass * T(1000)
    return -fac - (adt_slope > zero(T) ? adt_slope * scale : zero(T))
end

@inline function _enzyme_smooth_lambda_c3(
    fac,
    vcmax,
    stress,
    b,
    co2,
    temperature,
    apar,
    daylength,
    lpjmlparams,
    photoparams,
)
    T = typeof(fac)
    lambda = NeuralCrop.compute_lambda_c3_solution(
        fac,
        vcmax,
        stress,
        b,
        co2,
        temperature,
        apar,
        daylength,
        lpjmlparams,
        photoparams,
    )
    for _ in 1:8
        residual = fac * (one(T) - lambda) - NeuralCrop.c3_adtmm_scalar_impl(
            lambda,
            vcmax,
            stress,
            b,
            co2,
            temperature,
            apar,
            daylength,
            lpjmlparams,
            photoparams,
        )
        slope = _enzyme_c3_lambda_residual_slope(
            fac,
            lambda,
            vcmax,
            stress,
            b,
            co2,
            temperature,
            apar,
            daylength,
            lpjmlparams,
            photoparams,
        )
        lambda -= residual / slope
    end
    return lambda
end

@inline _enzyme_smooth_lambda(::Val{:C3}, args...) =
    _enzyme_smooth_lambda_c3(args...)

@inline function _enzyme_smooth_lambda(
    ::Val{:C4},
    fac::T,
    vcmax,
    stress,
    b,
    co2,
    temperature,
    apar,
    daylength,
    lpjmlparams,
    photoparams,
) where {T}
    lambda = NeuralCrop.compute_lambda_c4_solution(
        fac,
        vcmax,
        stress,
        b,
        temperature,
        apar,
        daylength,
        lpjmlparams,
        photoparams,
    )
    for _ in 1:8
        assimilation = NeuralCrop.c4_adtmm_scalar_impl(
            lambda,
            vcmax,
            stress,
            b,
            temperature,
            apar,
            daylength,
            lpjmlparams,
            photoparams,
        )
        phipi = min(one(T), lambda / T(photoparams.lambdamc4))
        light = stress * T(lpjmlparams.alphac4) * apar *
            T(photoparams.cmass) * T(photoparams.cq) / daylength
        je = phipi * light
        jc = vcmax / T(24)
        je_slope = lambda < T(photoparams.lambdamc4) ?
            light / T(photoparams.lambdamc4) : zero(T)
        total = je + jc
        root = sqrt(max(
            zero(T),
            total * total - T(4) * T(lpjmlparams.theta) * je * jc,
        ))
        root_slope = root > zero(T) ?
            (total - T(2) * T(lpjmlparams.theta) * jc) * je_slope / root :
            zero(T)
        gross_slope = (je_slope - root_slope) * daylength /
            (T(2) * T(lpjmlparams.theta))
        scale = (temperature + T(273.15)) / T(photoparams.p) * T(8.314) /
            T(photoparams.cmass) * T(1000)
        slope = -fac - (assimilation > zero(T) ? gross_slope * scale : zero(T))
        lambda -= (fac * (one(T) - lambda) - assimilation) / slope
    end
    return lambda
end

# The production lambda wrapper packs active CFT scalars into a NamedTuple
# before launching its kernel. Keep the AD path scalar and explicit so Enzyme
# can carry the CFT tangent through the continuous fixed-iteration solve.
function _enzyme_solve_lambda!(
    state::NeuralCrop.ModelState,
    cft::NeuralCrop.CFTParameters,
    pet,
    temperature,
    co2,
    lpjmlparams,
    photoparams,
    pathway::Union{Val{:C3}, Val{:C4}},
)
    T = eltype(NeuralCrop.crop_photosynthesis_auxiliary(state).lambda)
    lambda = NeuralCrop.crop_photosynthesis_auxiliary(state).lambda
    vcmax = NeuralCrop.crop_photosynthesis_auxiliary(state).vcmax
    temperature_stress = NeuralCrop.crop_photosynthesis_auxiliary(state).temperature_stress
    conductance = NeuralCrop.crop_canopy_auxiliary(state).canopy_conductance
    fpar = NeuralCrop.crop_canopy_auxiliary(state).fpar
    apar = NeuralCrop.crop_canopy_auxiliary(state).apar
    for cell in eachindex(lambda)
        co2_cell = co2[length(co2) == 1 ? 1 : cell]
        gpd, fac = NeuralCrop.compute_canopy_water_supply(
            pet.daylength[cell], conductance[cell], T(cft.gmin), fpar[cell], co2_cell,
        )
        lambda[cell] = if gpd > T(1e-5) &&
                          temperature_stress[cell] >= T(1e-2) &&
                          pet.daylength[cell] > zero(T) && co2_cell > zero(T)
            _enzyme_smooth_lambda(
                pathway,
                fac,
                vcmax[cell],
                temperature_stress[cell],
                T(cft.b),
                co2_cell,
                temperature[cell],
                apar[cell],
                pet.daylength[cell],
                lpjmlparams,
                photoparams,
            )
        else
            zero(T)
        end
    end
    return nothing
end

function _enzyme_evaporation!(state::NeuralCrop.ModelState, pet_eeq, lpjmlparams, layer_depth)
    T = eltype(pet_eeq)
    crop_fpar = NeuralCrop.crop_canopy_auxiliary(state).fpar
    crop_transpiration = NeuralCrop.crop_fluxes(state).water.transpiration_layer
    crop_canopy_wet = NeuralCrop.crop_canopy_auxiliary(state).canopy_wet
    soil_water = NeuralCrop.soil_water_auxiliary(state)
    soil_water_flux = NeuralCrop.soil_water_fluxes(state)
    surface_litter = NeuralCrop.soil_surface_litter_prognostic(state)
    surface_litter_aux = NeuralCrop.soil_surface_litter_auxiliary(state)
    surface_litter_flux = NeuralCrop.soil_surface_litter_fluxes(state)
    priestley_taylor = lpjmlparams.PRIESTLEY_TAYLOR
    for cell in eachindex(pet_eeq)
        evap_energy = pet_eeq[cell] * priestley_taylor *
            max(one(T) - crop_fpar[cell], T(0.05))
        crop_transpiration_sum = crop_transpiration[1, cell] +
            crop_transpiration[2, cell] + crop_transpiration[3, cell] +
            crop_transpiration[4, cell] + crop_transpiration[5, cell]
        available_evaporation = max(
            pet_eeq[cell] * priestley_taylor *
                (one(T) - crop_canopy_wet[cell]) - crop_transpiration_sum,
            zero(T),
        )
        surface_litter_flux.evaporation[cell] = NeuralCrop.compute_litter_evaporation(
            evap_energy,
            available_evaporation,
            surface_litter.water_storage[cell],
            surface_litter_aux.water_capacity[cell],
            surface_litter.cover[cell],
        )
        surface_litter.water_storage[cell] -= surface_litter_flux.evaporation[cell]

        liquid_1 = soil_water.relative_content[1, cell] *
            soil_water.holding_capacity_storage[1, cell] + soil_water.free_water[1, cell]
        liquid_2 = soil_water.relative_content[2, cell] *
            soil_water.holding_capacity_storage[2, cell] + soil_water.free_water[2, cell]
        liquid_3 = soil_water.relative_content[3, cell] *
            soil_water.holding_capacity_storage[3, cell] + soil_water.free_water[3, cell]
        liquid_4 = soil_water.relative_content[4, cell] *
            soil_water.holding_capacity_storage[4, cell] + soil_water.free_water[4, cell]
        liquid_5 = soil_water.relative_content[5, cell] *
            soil_water.holding_capacity_storage[5, cell] + soil_water.free_water[5, cell]
        evap_ratio = zero(T)
        if evap_energy > T(1e-5) && available_evaporation > T(1e-5)
            soil_depth_evap = T(lpjmlparams.soildepth_evap)
            fraction_1 = min(one(T), max(soil_depth_evap, zero(T)) / layer_depth[1])
            water_evap = max(liquid_1 - crop_transpiration[1, cell], zero(T)) * fraction_1
            holding_evap = soil_water.holding_capacity_storage[1, cell] * fraction_1
            soil_depth_evap -= layer_depth[1]
            fraction_2 = min(one(T), max(soil_depth_evap, zero(T)) / layer_depth[2])
            water_evap += max(liquid_2 - crop_transpiration[2, cell], zero(T)) * fraction_2
            holding_evap += soil_water.holding_capacity_storage[2, cell] * fraction_2
            soil_depth_evap -= layer_depth[2]
            fraction_3 = min(one(T), max(soil_depth_evap, zero(T)) / layer_depth[3])
            water_evap += max(liquid_3 - crop_transpiration[3, cell], zero(T)) * fraction_3
            holding_evap += soil_water.holding_capacity_storage[3, cell] * fraction_3
            soil_depth_evap -= layer_depth[3]
            fraction_4 = min(one(T), max(soil_depth_evap, zero(T)) / layer_depth[4])
            water_evap += max(liquid_4 - crop_transpiration[4, cell], zero(T)) * fraction_4
            holding_evap += soil_water.holding_capacity_storage[4, cell] * fraction_4
            soil_depth_evap -= layer_depth[4]
            fraction_5 = min(one(T), max(soil_depth_evap, zero(T)) / layer_depth[5])
            water_evap += max(liquid_5 - crop_transpiration[5, cell], zero(T)) * fraction_5
            holding_evap += soil_water.holding_capacity_storage[5, cell] * fraction_5
            evap_ratio = NeuralCrop.compute_soil_evaporation_ratio(
                evap_energy,
                available_evaporation,
                water_evap,
                holding_evap,
                surface_litter.cover[cell],
            )
        end

        soil_depth_evap = T(lpjmlparams.soildepth_evap)
        fraction_1 = min(one(T), max(soil_depth_evap, zero(T)) / layer_depth[1])
        soil_water_flux.evaporation[1, cell] = liquid_1 * evap_ratio * fraction_1
        soil_depth_evap -= layer_depth[1]
        fraction_2 = min(one(T), max(soil_depth_evap, zero(T)) / layer_depth[2])
        soil_water_flux.evaporation[2, cell] = liquid_2 * evap_ratio * fraction_2
        soil_depth_evap -= layer_depth[2]
        fraction_3 = min(one(T), max(soil_depth_evap, zero(T)) / layer_depth[3])
        soil_water_flux.evaporation[3, cell] = liquid_3 * evap_ratio * fraction_3
        soil_depth_evap -= layer_depth[3]
        fraction_4 = min(one(T), max(soil_depth_evap, zero(T)) / layer_depth[4])
        soil_water_flux.evaporation[4, cell] = liquid_4 * evap_ratio * fraction_4
        soil_depth_evap -= layer_depth[4]
        fraction_5 = min(one(T), max(soil_depth_evap, zero(T)) / layer_depth[5])
        soil_water_flux.evaporation[5, cell] = liquid_5 * evap_ratio * fraction_5
    end
    return nothing
end

function _enzyme_partition_soil_water_ice!(state::NeuralCrop.ModelState)
    water = NeuralCrop.soil_water_prognostic(state)
    auxiliary = NeuralCrop.soil_water_auxiliary(state)
    for cell in eachindex(NeuralCrop.crop_fluxes(state).carbon.gross_assimilation)
        _enzyme_partition_soil_water_ice_layer!(water, auxiliary, 1, cell)
        _enzyme_partition_soil_water_ice_layer!(water, auxiliary, 2, cell)
        _enzyme_partition_soil_water_ice_layer!(water, auxiliary, 3, cell)
        _enzyme_partition_soil_water_ice_layer!(water, auxiliary, 4, cell)
        _enzyme_partition_soil_water_ice_layer!(water, auxiliary, 5, cell)
    end
    return nothing
end

function _enzyme_partition_soil_water_ice_layer!(water, auxiliary, layer, cell)
    T = eltype(water.storage)
    total = max(water.storage[layer, cell] + water.ice_storage[layer, cell], zero(T))
    ice = clamp(water.ice_storage[layer, cell], zero(T), total)
    pwp_fraction, available_ice, free_ice, relative_water, free_water =
        NeuralCrop.lpjml_water_ice_partition(
            total,
            ice,
            auxiliary.wilting_storage[layer, cell],
            auxiliary.holding_capacity_storage[layer, cell],
        )
    water.ice_storage[layer, cell] = ice
    water.storage[layer, cell] = total - ice
    water.wilting_ice_fraction[layer, cell] = pwp_fraction
    water.available_ice_storage[layer, cell] = available_ice
    water.free_ice_storage[layer, cell] = free_ice
    auxiliary.relative_content[layer, cell] = relative_water
    auxiliary.free_water[layer, cell] = free_water
    return nothing
end

function _enzyme_soil_evapotranspiration!(
    state::NeuralCrop.ModelState;
    irrigation::Bool = false,
)
    water = NeuralCrop.soil_water_prognostic(state)
    storage = water.storage
    transpiration = NeuralCrop.crop_fluxes(state).water.transpiration_layer
    evaporation = NeuralCrop.soil_water_fluxes(state).evaporation
    if irrigation
        field_capacity = NeuralCrop.soil_water_auxiliary(state).field_capacity
        layer_depth = NeuralCrop.soil_properties(state).layer_depth
        ice_storage = water.ice_storage
        for cell in axes(storage, 2), layer in axes(storage, 1)
            target_storage = field_capacity[layer, cell] * layer_depth[layer]
            storage[layer, cell] = max(
                target_storage - ice_storage[layer, cell],
                zero(eltype(storage)),
            )
        end
    else
        for cell in eachindex(NeuralCrop.crop_fluxes(state).carbon.gross_assimilation)
            storage[1, cell] -= transpiration[1, cell] + evaporation[1, cell]
            storage[2, cell] -= transpiration[2, cell] + evaporation[2, cell]
            storage[3, cell] -= transpiration[3, cell] + evaporation[3, cell]
            storage[4, cell] -= transpiration[4, cell] + evaporation[4, cell]
            storage[5, cell] -= transpiration[5, cell] + evaporation[5, cell]
        end
    end
    _enzyme_partition_soil_water_ice!(state)
    return nothing
end

"""Keep the AD root profile linked to the active CFT beta_root parameter.

The production initializer evaluates this profile before the daily transition,
which is correct for ordinary runs but makes `beta_root` invisible when the
post-sowing state is treated as fixed. Reassigning the same profile inside the
opt-in AD path preserves baseline primal parity and exposes the intended
initial-condition sensitivity without changing production initialization.
"""
function _enzyme_apply_root_distribution!(state::NeuralCrop.ModelState, beta_root)
    distribution = state.inputs.crop.root.distribution
    T = eltype(distribution)
    beta = T(beta_root)
    beta20 = beta ^ 20
    beta50 = beta ^ 50
    beta100 = beta ^ 100
    beta200 = beta ^ 200
    beta300 = beta ^ 300
    total = one(T) - beta300
    distribution[1] = (one(T) - beta20) / total
    distribution[2] = (beta20 - beta50) / total
    distribution[3] = (beta50 - beta100) / total
    distribution[4] = (beta100 - beta200) / total
    distribution[5] = (beta200 - beta300) / total
    return nothing
end

"""Keep same-precision AD parameter bundles on the differentiable path."""
@inline function _enzyme_model_parameters(
    ::Type{T},
    parameters::NeuralCrop.ModelParameters{T},
) where {T <: AbstractFloat}
    return parameters
end

@inline function _enzyme_model_parameters(
    ::Type{T},
    parameters::NeuralCrop.ModelParameters,
) where {T <: AbstractFloat}
    return NeuralCrop.convert_precision(T, parameters)
end

@inline function _enzyme_pathway(cft::NeuralCrop.CFTParameters)
    cft.path == 1 && return Val(:C3)
    cft.path == 2 && return Val(:C4)
    throw(ArgumentError("unsupported crop photosynthesis path $(cft.path)"))
end

"""
    _enzyme_continuous_transition!(state, cft, global_parameters, climate, day, observable)

Run the existing daily process kernels between fixed discrete lifecycle events.
Sowing, harvest, failed-crop termination, and output recording are handled by
the outer station driver; this function is the opt-in continuous AD path.
"""
function _enzyme_continuous_transition!(
    state::NeuralCrop.ModelState,
    cft::NeuralCrop.CFTParameters,
    global_parameters,
    climate,
    day::Integer,
    observable::Symbol,
    layer_depth,
    irrigation::Bool,
    nitrogen_limit_vcmax::Bool,
    pathway::Union{Val{:C3}, Val{:C4}},
)
    T = eltype(NeuralCrop.crop_prognostic(state).canopy.lai)
    _enzyme_apply_root_distribution!(state, cft.beta_root)
    model_parameters = _enzyme_model_parameters(T, global_parameters)
    global_params = model_parameters.lpjml
    photo_params = model_parameters.photosynthesis
    snow_params = model_parameters.snow
    thermal_params = model_parameters.soil_thermal
    decomp_params = model_parameters.soil_decomposition
    climbuf = state.prognostic.climate
    pet = state.auxiliary.pet
    managed_land = state.inputs.management
    daily_weather = state.inputs.weather

    current_co2 = NeuralCrop.readclimate!(climate, daily_weather, day)
    NeuralCrop.update_climbuf!(
        cft,
        daily_weather.temp,
        climbuf,
        day;
        update_vernalization_requirement = false,
    )
    # The production driver applies its configured tillage treatment before
    # daily litter mixing. The fixed-event AD path uses the same daily
    # management treatment; sowing and harvest remain outer events.
    NeuralCrop.litter_tillage!(state, state)
    NeuralCrop.tillage_hydraulics!(state; lpjmlparams = global_params)
    NeuralCrop.litter_bioturbation!(state; lpjmlparams = global_params)

    NeuralCrop._pathway_albedo!(pathway, cft, state, state, pet, true)
    NeuralCrop.petpar!(
        pet,
        day % 365 == 0 ? 365 : day % 365,
        managed_land.latitude,
        daily_weather.temp,
        daily_weather.lwr,
        daily_weather.swr,
    )
    NeuralCrop.snow!(state, daily_weather; snowparams = snow_params, lpjmlparams = global_params)
    NeuralCrop.pedotransfer!(state; lpjmlparams = global_params)
    NeuralCrop.update_surface_litter_properties!(state; thermalparams = thermal_params)
    NeuralCrop.soil_temperature!(
        state,
        daily_weather.temp,
        climbuf.atemp_mean;
        thermalparams = thermal_params,
        snowparams = snow_params,
    )
    NeuralCrop.soil_cn_decomposition!(
        state;
        lpjmlparams = global_params,
        soil_decomp_params = decomp_params,
    )
    if hasproperty(climate, :no3_deposition) || hasproperty(climate, :nh4_deposition)
        NeuralCrop.nitrogen_deposition!(
            state,
            daily_weather.no3_deposition,
            daily_weather.nh4_deposition,
        )
    end
    NeuralCrop.phenology_crop!(
        state,
        climbuf.V_req,
        cft,
        daily_weather.temp,
        pet.daylength,
    )

    NeuralCrop.interception!(
        state,
        cft,
        pet.eeq,
        daily_weather.prec;
        lpjmlparams = global_params,
    )
    NeuralCrop.pedotransfer!(state; lpjmlparams = global_params)
    NeuralCrop.soil_infiltration!(
        state,
        state,
        daily_weather.prec;
        snowmelt = NeuralCrop.soil_snow_fluxes(state).melt,
        air_temperature = daily_weather.temp,
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
        cft,
        pet,
        state,
        daily_weather.temp;
        photoparams = photo_params,
    )
    NeuralCrop.photosynthesis!(
        pathway,
        cft,
        state,
        NeuralCrop.crop_canopy_auxiliary(state).apar,
        pet.daylength,
        daily_weather.temp,
        current_co2;
        comp_vcmax = true,
        lpjmlparams = global_params,
        photoparams = photo_params,
    )
    NeuralCrop.transpiration!(
        NeuralCrop.crop_fluxes(state).carbon.water_limited_assimilation,
        cft,
        state,
        pet,
        state,
        current_co2;
        lpjmlparams = global_params,
    )
    _enzyme_solve_lambda!(
        state,
        cft,
        pet,
        daily_weather.temp,
        current_co2,
        global_params,
        photo_params,
        pathway,
    )
    if nitrogen_limit_vcmax
        # Match the nitrogen-limited production path: derive demand and acquire
        # nitrogen from the potential capacity before applying the leaf-N Vcmax
        # constraint to the final photosynthesis evaluation.
        NeuralCrop.crop_nitrogen!(
            state,
            cft,
            state,
            NeuralCrop.crop_photosynthesis_auxiliary(state).potential_vcmax,
            daily_weather.temp;
            auto_fertilizer = false,
            lpjmlparams = global_params,
        )
        NeuralCrop.limit_vcmax_by_nitrogen!(
            state,
            cft,
            daily_weather.temp;
            lpjmlparams = global_params,
        )
    end
    NeuralCrop.photosynthesis!(
        pathway,
        cft,
        state,
        NeuralCrop.crop_canopy_auxiliary(state).apar,
        pet.daylength,
        daily_weather.temp,
        current_co2;
        comp_vcmax = false,
        lpjmlparams = global_params,
        photoparams = photo_params,
    )
    _enzyme_crop_carbon!(
        state,
        cft,
        daily_weather.temp,
        NeuralCrop.soil_thermal_prognostic(state).temperature,
        global_params,
    )
    if nitrogen_limit_vcmax
        # Carbon allocation changes organ weights; redistribute the acquired
        # total plant nitrogen exactly where the production driver does.
        NeuralCrop.allocate_crop_nitrogen!(state, cft)
    else
        NeuralCrop.crop_nitrogen!(
            state,
            cft,
            state,
            NeuralCrop.crop_photosynthesis_auxiliary(state).vcmax,
            daily_weather.temp;
            auto_fertilizer = false,
            lpjmlparams = global_params,
        )
    end
    _enzyme_evaporation!(state, pet.eeq, global_params, layer_depth)
    _enzyme_soil_evapotranspiration!(state; irrigation)
    NeuralCrop.post_crop_nitrogen_losses!(
        state;
        air_temperature = daily_weather.temp,
        wind_speed = daily_weather.wind,
        lpjmlparams = global_params,
    )
    return nothing
end

@generated function _replace_cft_parameters(
    base_cft::NeuralCrop.CFTParameters{T, S},
    theta::AbstractVector{T},
    ::Val{parameter_names},
) where {T, S, parameter_names}
    fields = fieldnames(NeuralCrop.CFTParameters)
    values = map(fields) do field
        if field === :temp_co2
            index = findfirst(==(:temp_co2_high), parameter_names)
            return index === nothing ?
                :(getfield(base_cft, :temp_co2)) :
                :(NeuralCrop.TempCO2{$T}(
                    getfield(getfield(base_cft, :temp_co2), :low),
                    theta[$index],
                ))
        elseif field === :temp_photos
            index = findfirst(==(:temp_photos_high), parameter_names)
            return index === nothing ?
                :(getfield(base_cft, :temp_photos)) :
                :(NeuralCrop.TempPhotos{$T}(
                    getfield(getfield(base_cft, :temp_photos), :low),
                    theta[$index],
                ))
        elseif field === :basetemp
            index = findfirst(==(:basetemp_low), parameter_names)
            return index === nothing ?
                :(getfield(base_cft, :basetemp)) :
                :(NeuralCrop.BaseTemp{$T}(
                    theta[$index],
                    getfield(getfield(base_cft, :basetemp), :high),
                ))
        end
        index = findfirst(==(field), parameter_names)
        index === nothing ? :(getfield(base_cft, $(QuoteNode(field)))) : :(theta[$index])
    end
    return :(NeuralCrop.CFTParameters{$T, $S}($(values...)))
end

@generated function _replace_cft_shadow(
    zero_cft::NeuralCrop.CFTParameters{T, S},
    dtheta::AbstractVector{T},
    ::Val{parameter_names},
) where {T, S, parameter_names}
    fields = fieldnames(NeuralCrop.CFTParameters)
    values = map(fields) do field
        if field === :temp_co2
            index = findfirst(==(:temp_co2_high), parameter_names)
            return index === nothing ?
                :(getfield(zero_cft, :temp_co2)) :
                :(NeuralCrop.TempCO2{$T}(
                    getfield(getfield(zero_cft, :temp_co2), :low),
                    dtheta[$index],
                ))
        elseif field === :temp_photos
            index = findfirst(==(:temp_photos_high), parameter_names)
            return index === nothing ?
                :(getfield(zero_cft, :temp_photos)) :
                :(NeuralCrop.TempPhotos{$T}(
                    getfield(getfield(zero_cft, :temp_photos), :low),
                    dtheta[$index],
                ))
        elseif field === :basetemp
            index = findfirst(==(:basetemp_low), parameter_names)
            return index === nothing ?
                :(getfield(zero_cft, :basetemp)) :
                :(NeuralCrop.BaseTemp{$T}(
                    dtheta[$index],
                    getfield(getfield(zero_cft, :basetemp), :high),
                ))
        end
        index = findfirst(==(field), parameter_names)
        index === nothing ? :(getfield(zero_cft, $(QuoteNode(field)))) : :(dtheta[$index])
    end
    return :(NeuralCrop.CFTParameters{$T, $S}($(values...)))
end

_replace_cft_shadow(zero_cft, dtheta, parameter_names::Tuple) =
    _replace_cft_shadow(zero_cft, dtheta, Val(parameter_names))

function Enzyme.EnzymeRules.forward(
    config::Enzyme.EnzymeRules.FwdConfig,
    ::Enzyme.Const{typeof(_replace_cft_parameters)},
    ::Type,
    base_cft::Enzyme.Const,
    theta::Enzyme.Duplicated,
    parameter_names::Enzyme.Const,
)
    primal = _replace_cft_parameters(base_cft.val, theta.val, parameter_names.val)
    shadow = _replace_cft_shadow(
        Enzyme.make_zero(base_cft.val),
        theta.dval,
        parameter_names.val,
    )
    if Enzyme.EnzymeRules.needs_primal(config) && Enzyme.EnzymeRules.needs_shadow(config)
        return Enzyme.Duplicated(primal, shadow)
    elseif Enzyme.EnzymeRules.needs_shadow(config)
        return shadow
    elseif Enzyme.EnzymeRules.needs_primal(config)
        return primal
    end
    return nothing
end

function _replace_cft_parameters(
    base_cft::NeuralCrop.CFTParameters{T, S},
    theta::AbstractVector{T},
    parameter_names::Tuple,
)::NeuralCrop.CFTParameters{T, S} where {T, S}
    length(theta) == length(parameter_names) || throw(DimensionMismatch(
        "theta and parameter_names must have identical lengths",
    ))
    return _replace_cft_parameters(base_cft, theta, Val(parameter_names))
end

"""Replace only the selected soil-related fields in `LPJmLParams`."""
@generated function _replace_lpjml_parameters(
    base_lpjml::NeuralCrop.LPJmLParams{T},
    theta::AbstractVector{T},
    ::Val{parameter_names},
) where {T, parameter_names}
    fields = fieldnames(NeuralCrop.LPJmLParams)
    values = map(fields) do field
        if field === :k_soil10
            index = findfirst(==(:k_soil10_fast), parameter_names)
            fast = index === nothing ?
                :(getfield(getfield(base_lpjml, :k_soil10), :fast)) :
                :(theta[$index])
            return :(NeuralCrop.K_Soil10{$T}(
                $fast,
                getfield(getfield(base_lpjml, :k_soil10), :slow),
            ))
        end
        index = findfirst(==(field), parameter_names)
        return index === nothing ?
            :(getfield(base_lpjml, $(QuoteNode(field)))) :
            :(theta[$index])
    end
    return :(NeuralCrop.LPJmLParams{$T}($(values...)))
end

@generated function _replace_lpjml_shadow(
    zero_lpjml::NeuralCrop.LPJmLParams{T},
    dtheta::AbstractVector{T},
    ::Val{parameter_names},
) where {T, parameter_names}
    fields = fieldnames(NeuralCrop.LPJmLParams)
    values = map(fields) do field
        if field === :k_soil10
            index = findfirst(==(:k_soil10_fast), parameter_names)
            fast = index === nothing ?
                :(getfield(getfield(zero_lpjml, :k_soil10), :fast)) :
                :(dtheta[$index])
            return :(NeuralCrop.K_Soil10{$T}(
                $fast,
                getfield(getfield(zero_lpjml, :k_soil10), :slow),
            ))
        end
        index = findfirst(==(field), parameter_names)
        return index === nothing ?
            :(getfield(zero_lpjml, $(QuoteNode(field)))) :
            :(dtheta[$index])
    end
    return :(NeuralCrop.LPJmLParams{$T}($(values...)))
end

_replace_lpjml_shadow(zero_lpjml, dtheta, parameter_names::Tuple) =
    _replace_lpjml_shadow(zero_lpjml, dtheta, Val(parameter_names))

function Enzyme.EnzymeRules.forward(
    config::Enzyme.EnzymeRules.FwdConfig,
    ::Enzyme.Const{typeof(_replace_lpjml_parameters)},
    ::Type,
    base_lpjml::Enzyme.Const,
    theta::Enzyme.Duplicated,
    parameter_names::Enzyme.Const,
)
    primal = _replace_lpjml_parameters(
        base_lpjml.val, theta.val, parameter_names.val,
    )
    shadow = _replace_lpjml_shadow(
        Enzyme.make_zero(base_lpjml.val),
        theta.dval,
        parameter_names.val,
    )
    if Enzyme.EnzymeRules.needs_primal(config) && Enzyme.EnzymeRules.needs_shadow(config)
        return Enzyme.Duplicated(primal, shadow)
    elseif Enzyme.EnzymeRules.needs_shadow(config)
        return shadow
    elseif Enzyme.EnzymeRules.needs_primal(config)
        return primal
    end
    return nothing
end

function _replace_lpjml_parameters(
    base_lpjml::NeuralCrop.LPJmLParams{T},
    theta::AbstractVector{T},
    parameter_names::Tuple,
)::NeuralCrop.LPJmLParams{T} where {T}
    length(theta) == length(parameter_names) || throw(DimensionMismatch(
        "theta and parameter_names must have identical lengths",
    ))
    return _replace_lpjml_parameters(base_lpjml, theta, Val(parameter_names))
end

"""Rebuild ModelParameters while keeping all non-soil bundles fixed."""
@generated function _replace_model_parameters(
    base_parameters::NeuralCrop.ModelParameters{T},
    theta::AbstractVector{T},
    ::Val{parameter_names},
) where {T, parameter_names}
    return :(NeuralCrop.ModelParameters{$T}(
        _replace_lpjml_parameters(
            getfield(base_parameters, :lpjml),
            theta,
            Val($(QuoteNode(parameter_names))),
        ),
        getfield(base_parameters, :photosynthesis),
        getfield(base_parameters, :snow),
        getfield(base_parameters, :soil_thermal),
        getfield(base_parameters, :soil_decomposition),
    ))
end

@generated function _replace_model_shadow(
    zero_parameters::NeuralCrop.ModelParameters{T},
    dtheta::AbstractVector{T},
    ::Val{parameter_names},
) where {T, parameter_names}
    return :(NeuralCrop.ModelParameters{$T}(
        _replace_lpjml_shadow(
            getfield(zero_parameters, :lpjml),
            dtheta,
            Val($(QuoteNode(parameter_names))),
        ),
        getfield(zero_parameters, :photosynthesis),
        getfield(zero_parameters, :snow),
        getfield(zero_parameters, :soil_thermal),
        getfield(zero_parameters, :soil_decomposition),
    ))
end

_replace_model_shadow(zero_parameters, dtheta, parameter_names::Tuple) =
    _replace_model_shadow(zero_parameters, dtheta, Val(parameter_names))

function Enzyme.EnzymeRules.forward(
    config::Enzyme.EnzymeRules.FwdConfig,
    ::Enzyme.Const{typeof(_replace_model_parameters)},
    ::Type,
    base_parameters::Enzyme.Const,
    theta::Enzyme.Duplicated,
    parameter_names::Enzyme.Const,
)
    primal = _replace_model_parameters(
        base_parameters.val, theta.val, parameter_names.val,
    )
    shadow = _replace_model_shadow(
        Enzyme.make_zero(base_parameters.val),
        theta.dval,
        parameter_names.val,
    )
    if Enzyme.EnzymeRules.needs_primal(config) && Enzyme.EnzymeRules.needs_shadow(config)
        return Enzyme.Duplicated(primal, shadow)
    elseif Enzyme.EnzymeRules.needs_shadow(config)
        return shadow
    elseif Enzyme.EnzymeRules.needs_primal(config)
        return primal
    end
    return nothing
end

function _replace_model_parameters(
    base_parameters::NeuralCrop.ModelParameters{T},
    theta::AbstractVector{T},
    parameter_names::Tuple,
)::NeuralCrop.ModelParameters{T} where {T}
    length(theta) == length(parameter_names) || throw(DimensionMismatch(
        "theta and parameter_names must have identical lengths",
    ))
    return _replace_model_parameters(base_parameters, theta, Val(parameter_names))
end

function _daily_transition_observable(state::NeuralCrop.ModelState, observable::Symbol)
    crop_flux = NeuralCrop.crop_fluxes(state)
    soil_carbon_flux = NeuralCrop.soil_carbon_fluxes(state)
    soil_water_flux = NeuralCrop.soil_water_fluxes(state)
    crop_water_flux = crop_flux.water
    litter_flux = NeuralCrop.soil_surface_litter_fluxes(state)
    if observable === :gpp
        return crop_flux.carbon.gross_assimilation[1]
    elseif observable === :reco
        return crop_flux.carbon.respiration[1] +
            crop_flux.carbon.leaf_respiration[1] +
            soil_carbon_flux.heterotrophic_respiration[1]
    elseif observable === :et
        value = crop_water_flux.interception[1] + litter_flux.evaporation[1]
        value += crop_water_flux.transpiration_layer[1, 1] + soil_water_flux.evaporation[1, 1]
        value += crop_water_flux.transpiration_layer[2, 1] + soil_water_flux.evaporation[2, 1]
        value += crop_water_flux.transpiration_layer[3, 1] + soil_water_flux.evaporation[3, 1]
        value += crop_water_flux.transpiration_layer[4, 1] + soil_water_flux.evaporation[4, 1]
        value += crop_water_flux.transpiration_layer[5, 1] + soil_water_flux.evaporation[5, 1]
        return value
    end
    throw(ArgumentError("unsupported daily transition observable $observable"))
end

"""
    NeuralCrop.enzyme_seasonal_loss(
        theta, state, base_cft, global_parameters, climate,
        parameter_names, days, layer_depth, context,
    )

Advance one fixed daily process trajectory and accumulate normalized mean
squared losses for GPP, RECO, and ET. The context masks only decide which
observations contribute; every day still advances the full continuous state.
Discrete events remain outside this objective.
"""
function NeuralCrop.enzyme_seasonal_loss(
    theta::AbstractVector{T},
    state::NeuralCrop.ModelState,
    base_cft::NeuralCrop.CFTParameters,
    global_parameters,
    climate,
    parameter_names::Tuple,
    days::Tuple,
    layer_depth,
    context::NeuralCrop.ADSeasonContext;
    irrigation::Bool = false,
) where {T <: AbstractFloat}
    pathway = _enzyme_pathway(base_cft)
    cft = _replace_cft_parameters(base_cft, theta, parameter_names)
    gpp_loss = zero(T)
    reco_loss = zero(T)
    et_loss = zero(T)
    day_index = 1
    while day_index <= length(days)
        _enzyme_continuous_transition!(
            state,
            cft,
            global_parameters,
            climate,
            days[day_index],
            :et,
            layer_depth,
            irrigation,
            false,
            pathway,
        )
        if context.growth_mask[day_index]
            crop_flux = NeuralCrop.crop_fluxes(state)
            soil_carbon_flux = NeuralCrop.soil_carbon_fluxes(state)
            crop_water_flux = crop_flux.water
            soil_water_flux = NeuralCrop.soil_water_fluxes(state)
            litter_flux = NeuralCrop.soil_surface_litter_fluxes(state)
            gpp = crop_flux.carbon.gross_assimilation[1]
            reco = crop_flux.carbon.respiration[1] +
                crop_flux.carbon.leaf_respiration[1] +
                soil_carbon_flux.heterotrophic_respiration[1]
            et = crop_water_flux.interception[1] + litter_flux.evaporation[1]
            et += crop_water_flux.transpiration_layer[1, 1] + soil_water_flux.evaporation[1, 1]
            et += crop_water_flux.transpiration_layer[2, 1] + soil_water_flux.evaporation[2, 1]
            et += crop_water_flux.transpiration_layer[3, 1] + soil_water_flux.evaporation[3, 1]
            et += crop_water_flux.transpiration_layer[4, 1] + soil_water_flux.evaporation[4, 1]
            et += crop_water_flux.transpiration_layer[5, 1] + soil_water_flux.evaporation[5, 1]

            if context.valid_masks.gpp[day_index]
                residual = (gpp - T(context.observations.gpp[day_index])) /
                    T(context.scales.gpp)
                gpp_loss += residual * residual
            end
            if context.valid_masks.reco[day_index]
                residual = (reco - T(context.observations.reco[day_index])) /
                    T(context.scales.reco)
                reco_loss += residual * residual
            end
            if context.valid_masks.et[day_index]
                residual = (et - T(context.observations.et[day_index])) /
                    T(context.scales.et)
                et_loss += residual * residual
            end
        end
        day_index += 1
    end
    gpp_count = T(max(context.counts.gpp, 1))
    reco_count = T(max(context.counts.reco, 1))
    et_count = T(max(context.counts.et, 1))
    return gpp_loss / gpp_count + reco_loss / reco_count + et_loss / et_count
end

mutable struct _ADSeasonalLossAccumulator{T}
    gpp::T
    reco::T
    et::T
end

@inline function _enzyme_add_seasonal_loss!(
    losses::_ADSeasonalLossAccumulator{T},
    state::NeuralCrop.ModelState,
    index::Int,
    context::NeuralCrop.ADSeasonContext,
) where {T <: AbstractFloat}
    context.growth_mask[index] || return nothing

    crop_flux = NeuralCrop.crop_fluxes(state)
    soil_carbon_flux = NeuralCrop.soil_carbon_fluxes(state)
    crop_water_flux = crop_flux.water
    soil_water_flux = NeuralCrop.soil_water_fluxes(state)
    litter_flux = NeuralCrop.soil_surface_litter_fluxes(state)
    gpp = crop_flux.carbon.gross_assimilation[1]
    reco = crop_flux.carbon.respiration[1] +
        crop_flux.carbon.leaf_respiration[1] +
        soil_carbon_flux.heterotrophic_respiration[1]
    et = crop_water_flux.interception[1] + litter_flux.evaporation[1]
    et += crop_water_flux.transpiration_layer[1, 1] + soil_water_flux.evaporation[1, 1]
    et += crop_water_flux.transpiration_layer[2, 1] + soil_water_flux.evaporation[2, 1]
    et += crop_water_flux.transpiration_layer[3, 1] + soil_water_flux.evaporation[3, 1]
    et += crop_water_flux.transpiration_layer[4, 1] + soil_water_flux.evaporation[4, 1]
    et += crop_water_flux.transpiration_layer[5, 1] + soil_water_flux.evaporation[5, 1]

    if context.valid_masks.gpp[index]
        residual = (gpp - T(context.observations.gpp[index])) / T(context.scales.gpp)
        losses.gpp += residual * residual
    end
    if context.valid_masks.reco[index]
        residual = (reco - T(context.observations.reco[index])) / T(context.scales.reco)
        losses.reco += residual * residual
    end
    if context.valid_masks.et[index]
        residual = (et - T(context.observations.et[index])) / T(context.scales.et)
        losses.et += residual * residual
    end
    return nothing
end

@inline function _enzyme_seasonal_loss_value(
    losses::_ADSeasonalLossAccumulator{T},
    context::NeuralCrop.ADSeasonContext,
) where {T <: AbstractFloat}
    gpp_count = T(max(context.counts.gpp, 1))
    reco_count = T(max(context.counts.reco, 1))
    et_count = T(max(context.counts.et, 1))
    return losses.gpp / gpp_count + losses.reco / reco_count + losses.et / et_count
end

"""Run and differentiate one fixed seasonal block.

The context counts remain those of the full season, so summing block losses
reproduces the normalization used by `enzyme_seasonal_loss`.
"""
function _enzyme_seasonal_loss_block(
    theta::AbstractVector{T},
    state::NeuralCrop.ModelState,
    base_cft::NeuralCrop.CFTParameters,
    global_parameters,
    climate,
    parameter_names::Tuple,
    days::Tuple,
    layer_depth,
    context::NeuralCrop.ADSeasonContext,
    day_range::UnitRange{Int},
    irrigation::Bool,
    pathway::Union{Val{:C3}, Val{:C4}},
) where {T <: AbstractFloat}
    cft = _replace_cft_parameters(base_cft, theta, parameter_names)
    losses = _ADSeasonalLossAccumulator(zero(T), zero(T), zero(T))
    for index in day_range
        _enzyme_continuous_transition!(
            state,
            cft,
            global_parameters,
            climate,
            days[index],
            :et,
            layer_depth,
            irrigation,
            false,
            pathway,
        )
        _enzyme_add_seasonal_loss!(losses, state, index, context)
    end
    return _enzyme_seasonal_loss_value(losses, context)
end

"""
    NeuralCrop.enzyme_seasonal_soil_loss(
        theta_soil, state, base_cft, base_parameters, soil_parameter_names,
        climate, days, layer_depth, context,
    )

Accumulate the fixed-event seasonal loss while differentiating only selected
soil fields inside `ModelParameters`. The CFT and all non-soil parameter
bundles remain constant.
"""
function NeuralCrop.enzyme_seasonal_soil_loss(
    theta_soil::AbstractVector{T},
    state::NeuralCrop.ModelState,
    base_cft::NeuralCrop.CFTParameters,
    base_parameters::NeuralCrop.ModelParameters,
    soil_parameter_names::Tuple,
    climate,
    days::Tuple,
    layer_depth,
    context::NeuralCrop.ADSeasonContext;
    irrigation::Bool = false,
) where {T <: AbstractFloat}
    pathway = _enzyme_pathway(base_cft)
    model_parameters = _replace_model_parameters(
        base_parameters, theta_soil, soil_parameter_names,
    )
    losses = _ADSeasonalLossAccumulator(zero(T), zero(T), zero(T))
    for index in eachindex(days)
        _enzyme_continuous_transition!(
            state,
            base_cft,
            model_parameters,
            climate,
            days[index],
            :et,
            layer_depth,
            irrigation,
            false,
            pathway,
        )
        _enzyme_add_seasonal_loss!(losses, state, index, context)
    end
    return _enzyme_seasonal_loss_value(losses, context)
end

"""
    NeuralCrop.enzyme_seasonal_joint_loss(
        theta_cft, theta_soil, state, base_cft, base_parameters,
        cft_parameter_names, soil_parameter_names, climate, days,
        layer_depth, context,
    )

Accumulate one fixed-event seasonal loss while differentiating selected CFT
and soil fields together. Existing CFT-only and soil-only entry points remain
unchanged.
"""
function NeuralCrop.enzyme_seasonal_joint_loss(
    theta_cft::AbstractVector{T},
    theta_soil::AbstractVector{T},
    state::NeuralCrop.ModelState,
    base_cft::NeuralCrop.CFTParameters,
    base_parameters::NeuralCrop.ModelParameters,
    cft_parameter_names::Tuple,
    soil_parameter_names::Tuple,
    climate,
    days::Tuple,
    layer_depth,
    context::NeuralCrop.ADSeasonContext;
    irrigation::Bool = false,
) where {T <: AbstractFloat}
    pathway = _enzyme_pathway(base_cft)
    cft = _replace_cft_parameters(base_cft, theta_cft, cft_parameter_names)
    model_parameters = _replace_model_parameters(
        base_parameters, theta_soil, soil_parameter_names,
    )
    losses = _ADSeasonalLossAccumulator(zero(T), zero(T), zero(T))
    for index in eachindex(days)
        _enzyme_continuous_transition!(
            state,
            cft,
            model_parameters,
            climate,
            days[index],
            :et,
            layer_depth,
            irrigation,
            false,
            pathway,
        )
        _enzyme_add_seasonal_loss!(losses, state, index, context)
    end
    return _enzyme_seasonal_loss_value(losses, context)
end

function _enzyme_seasonal_joint_loss_block(
    theta_cft::AbstractVector{T},
    theta_soil::AbstractVector{T},
    state::NeuralCrop.ModelState,
    base_cft::NeuralCrop.CFTParameters,
    base_parameters::NeuralCrop.ModelParameters,
    cft_parameter_names::Tuple,
    soil_parameter_names::Tuple,
    climate,
    days::Tuple,
    layer_depth,
    context::NeuralCrop.ADSeasonContext,
    day_range::UnitRange{Int},
    irrigation::Bool,
    pathway::Union{Val{:C3}, Val{:C4}},
) where {T <: AbstractFloat}
    cft = _replace_cft_parameters(base_cft, theta_cft, cft_parameter_names)
    model_parameters = _replace_model_parameters(
        base_parameters, theta_soil, soil_parameter_names,
    )
    losses = _ADSeasonalLossAccumulator(zero(T), zero(T), zero(T))
    for index in day_range
        _enzyme_continuous_transition!(
            state,
            cft,
            model_parameters,
            climate,
            days[index],
            :et,
            layer_depth,
            irrigation,
            false,
            pathway,
        )
        _enzyme_add_seasonal_loss!(losses, state, index, context)
    end
    return _enzyme_seasonal_loss_value(losses, context)
end

function _enzyme_seasonal_soil_loss_block(
    theta_soil::AbstractVector{T},
    state::NeuralCrop.ModelState,
    base_cft::NeuralCrop.CFTParameters,
    base_parameters::NeuralCrop.ModelParameters,
    soil_parameter_names::Tuple,
    climate,
    days::Tuple,
    layer_depth,
    context::NeuralCrop.ADSeasonContext,
    day_range::UnitRange{Int},
    irrigation::Bool,
    pathway::Union{Val{:C3}, Val{:C4}},
) where {T <: AbstractFloat}
    model_parameters = _replace_model_parameters(
        base_parameters, theta_soil, soil_parameter_names,
    )
    losses = _ADSeasonalLossAccumulator(zero(T), zero(T), zero(T))
    for index in day_range
        _enzyme_continuous_transition!(
            state,
            base_cft,
            model_parameters,
            climate,
            days[index],
            :et,
            layer_depth,
            irrigation,
            false,
            pathway,
        )
        _enzyme_add_seasonal_loss!(losses, state, index, context)
    end
    return _enzyme_seasonal_loss_value(losses, context)
end

"""Compute a checkpointed reverse gradient for the independent soil block."""
function NeuralCrop.enzyme_seasonal_soil_gradient_blockwise(
    theta_soil::AbstractVector{T},
    state_factory::F,
    base_cft::NeuralCrop.CFTParameters,
    base_parameters::NeuralCrop.ModelParameters,
    soil_parameter_names::Tuple,
    climate,
    days::Tuple,
    layer_depth,
    context::NeuralCrop.ADSeasonContext;
    block_days::Integer = 30,
    irrigation::Bool = false,
) where {T <: AbstractFloat, F}
    length(days) == length(context.growth_mask) || throw(DimensionMismatch(
        "days and context must have identical lengths",
    ))
    ranges = _seasonal_block_ranges(days, Int(block_days))
    isempty(ranges) && throw(ArgumentError("days must not be empty"))
    pathway = _enzyme_pathway(base_cft)

    state, _ = _state_and_shadow(state_factory)
    NeuralCrop.enzyme_prepare_daily_state!(state)
    snapshots = Vector{typeof(state)}(undef, length(ranges))
    forward_primal = zero(T)
    for (block_index, day_range) in enumerate(ranges)
        snapshots[block_index] = deepcopy(state)
        forward_primal += _enzyme_seasonal_soil_loss_block(
            theta_soil,
            state,
            base_cft,
            base_parameters,
            soil_parameter_names,
            climate,
            days,
            layer_depth,
            context,
            day_range,
            irrigation,
            pathway,
        )
    end

    gradient = zeros(T, length(theta_soil))
    state_cotangent = nothing
    reverse_primal = zero(T)
    for block_index in length(ranges):-1:1
        block_state = deepcopy(snapshots[block_index])
        block_shadow = state_cotangent === nothing ?
            NeuralCrop.enzyme_zero_tangent(block_state) : deepcopy(state_cotangent)
        block_gradient = zeros(T, length(theta_soil))
        day_range = ranges[block_index]
        result = Enzyme.autodiff(
            Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal),
            _enzyme_seasonal_soil_loss_block,
            Enzyme.Duplicated(theta_soil, block_gradient),
            Enzyme.Duplicated(block_state, block_shadow),
            Enzyme.Const(base_cft),
            Enzyme.Const(base_parameters),
            Enzyme.Const(soil_parameter_names),
            Enzyme.Const(climate),
            Enzyme.Const(days),
            Enzyme.Const(layer_depth),
            Enzyme.Const(context),
            Enzyme.Const(day_range),
            Enzyme.Const(irrigation),
            Enzyme.Const(pathway),
        )
        reverse_primal += result[2]
        gradient .+= block_gradient
        state_cotangent = block_shadow
    end

    return (
        primal = reverse_primal,
        forward_primal = forward_primal,
        gradient,
        block_days = Int(block_days),
        block_ranges = ranges,
    )
end

"""Compute one checkpointed reverse gradient for CFT and soil parameters."""
function NeuralCrop.enzyme_seasonal_joint_gradient_blockwise(
    theta_cft::AbstractVector{T},
    theta_soil::AbstractVector{T},
    state_factory::F,
    base_cft::NeuralCrop.CFTParameters,
    base_parameters::NeuralCrop.ModelParameters,
    cft_parameter_names::Tuple,
    soil_parameter_names::Tuple,
    climate,
    days::Tuple,
    layer_depth,
    context::NeuralCrop.ADSeasonContext;
    block_days::Integer = 30,
    irrigation::Bool = false,
) where {T <: AbstractFloat, F}
    length(theta_cft) == length(cft_parameter_names) || throw(DimensionMismatch(
        "CFT parameter values and names must have identical lengths",
    ))
    length(theta_soil) == length(soil_parameter_names) || throw(DimensionMismatch(
        "soil parameter values and names must have identical lengths",
    ))
    length(days) == length(context.growth_mask) || throw(DimensionMismatch(
        "days and context must have identical lengths",
    ))
    ranges = _seasonal_block_ranges(days, Int(block_days))
    isempty(ranges) && throw(ArgumentError("days must not be empty"))
    pathway = _enzyme_pathway(base_cft)

    state, _ = _state_and_shadow(state_factory)
    NeuralCrop.enzyme_prepare_daily_state!(state)
    snapshots = Vector{typeof(state)}(undef, length(ranges))
    forward_primal = zero(T)
    for (block_index, day_range) in enumerate(ranges)
        snapshots[block_index] = deepcopy(state)
        forward_primal += _enzyme_seasonal_joint_loss_block(
            theta_cft,
            theta_soil,
            state,
            base_cft,
            base_parameters,
            cft_parameter_names,
            soil_parameter_names,
            climate,
            days,
            layer_depth,
            context,
            day_range,
            irrigation,
            pathway,
        )
    end

    cft_gradient = zeros(T, length(theta_cft))
    soil_gradient = zeros(T, length(theta_soil))
    state_cotangent = nothing
    reverse_primal = zero(T)
    for block_index in length(ranges):-1:1
        block_state = deepcopy(snapshots[block_index])
        block_shadow = state_cotangent === nothing ?
            NeuralCrop.enzyme_zero_tangent(block_state) : deepcopy(state_cotangent)
        block_cft_gradient = zeros(T, length(theta_cft))
        block_soil_gradient = zeros(T, length(theta_soil))
        day_range = ranges[block_index]
        result = Enzyme.autodiff(
            Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal),
            _enzyme_seasonal_joint_loss_block,
            Enzyme.Duplicated(theta_cft, block_cft_gradient),
            Enzyme.Duplicated(theta_soil, block_soil_gradient),
            Enzyme.Duplicated(block_state, block_shadow),
            Enzyme.Const(base_cft),
            Enzyme.Const(base_parameters),
            Enzyme.Const(cft_parameter_names),
            Enzyme.Const(soil_parameter_names),
            Enzyme.Const(climate),
            Enzyme.Const(days),
            Enzyme.Const(layer_depth),
            Enzyme.Const(context),
            Enzyme.Const(day_range),
            Enzyme.Const(irrigation),
            Enzyme.Const(pathway),
        )
        reverse_primal += result[2]
        cft_gradient .+= block_cft_gradient
        soil_gradient .+= block_soil_gradient
        state_cotangent = block_shadow
    end

    return (;
        primal = reverse_primal,
        forward_primal,
        cft_gradient,
        soil_gradient,
        block_days = Int(block_days),
        block_ranges = ranges,
    )
end

function _seasonal_block_ranges(days, block_days::Int)
    block_days > 0 || throw(ArgumentError("block_days must be positive"))
    ranges = UnitRange{Int}[]
    first_day = 1
    while first_day <= length(days)
        last_day = min(first_day + block_days - 1, length(days))
        push!(ranges, first_day:last_day)
        first_day = last_day + 1
    end
    return ranges
end

"""
    NeuralCrop.enzyme_seasonal_gradient_blockwise(
        theta, state_factory, base_cft, global_parameters, climate,
        parameter_names, days, layer_depth, context; block_days=30,
    )

Compute the exact reverse gradient of a fixed seasonal loss with explicit
periodic state checkpoints. A forward pass stores one state snapshot at each
block boundary. Blocks are then differentiated from the final block backwards;
the state shadow from a later block seeds the output-state cotangent of the
preceding block. This is intentionally separate from the production runtime.
"""
function NeuralCrop.enzyme_seasonal_gradient_blockwise(
    theta::AbstractVector{T},
    state_factory::F,
    base_cft::NeuralCrop.CFTParameters,
    global_parameters,
    climate,
    parameter_names::Tuple,
    days::Tuple,
    layer_depth,
    context::NeuralCrop.ADSeasonContext;
    block_days::Integer = 30,
    irrigation::Bool = false,
) where {T <: AbstractFloat, F}
    length(days) == length(context.growth_mask) || throw(DimensionMismatch(
        "days and context must have identical lengths",
    ))
    ranges = _seasonal_block_ranges(days, Int(block_days))
    isempty(ranges) && throw(ArgumentError("days must not be empty"))
    pathway = _enzyme_pathway(base_cft)

    state, _ = _state_and_shadow(state_factory)
    NeuralCrop.enzyme_prepare_daily_state!(state)
    snapshots = Vector{typeof(state)}(undef, length(ranges))
    forward_primal = zero(T)
    for (block_index, day_range) in enumerate(ranges)
        snapshots[block_index] = deepcopy(state)
        forward_primal += _enzyme_seasonal_loss_block(
            theta,
            state,
            base_cft,
            global_parameters,
            climate,
            parameter_names,
            days,
            layer_depth,
            context,
            day_range,
            irrigation,
            pathway,
        )
    end

    gradient = zeros(T, length(theta))
    state_cotangent = nothing
    reverse_primal = zero(T)
    for block_index in length(ranges):-1:1
        block_state = deepcopy(snapshots[block_index])
        block_shadow = state_cotangent === nothing ?
            NeuralCrop.enzyme_zero_tangent(block_state) : deepcopy(state_cotangent)
        block_gradient = zeros(T, length(theta))
        day_range = ranges[block_index]
        result = Enzyme.autodiff(
            Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal),
            _enzyme_seasonal_loss_block,
            Enzyme.Duplicated(theta, block_gradient),
            Enzyme.Duplicated(block_state, block_shadow),
            Enzyme.Const(base_cft),
            Enzyme.Const(global_parameters),
            Enzyme.Const(climate),
            Enzyme.Const(parameter_names),
            Enzyme.Const(days),
            Enzyme.Const(layer_depth),
            Enzyme.Const(context),
            Enzyme.Const(day_range),
            Enzyme.Const(irrigation),
            Enzyme.Const(pathway),
        )
        reverse_primal += result[2]
        gradient .+= block_gradient
        state_cotangent = block_shadow
    end

    return (
        primal = reverse_primal,
        forward_primal = forward_primal,
        gradient = gradient,
        block_days = Int(block_days),
        block_ranges = ranges,
    )
end

"""
    NeuralCrop.enzyme_daily_transition_objective(
        theta, state, base_cft, global_parameters, climate,
        parameter_names, day, observable, layer_depth,
    )

Run the opt-in continuous portion of one existing C3 daily transition and
return one scalar daily observable. The caller should prepare the lifecycle
state with discrete events already fixed; output diagnostics are intentionally
excluded from this AD path.
"""
function NeuralCrop.enzyme_daily_transition_objective(
    theta::AbstractVector{T},
    state::NeuralCrop.ModelState,
    base_cft::NeuralCrop.CFTParameters,
    global_parameters,
    climate,
    parameter_names::Tuple,
    day::Integer,
    observable::Symbol,
    layer_depth;
    irrigation::Bool = false,
) where {T <: AbstractFloat}
    pathway = _enzyme_pathway(base_cft)
    cft = _replace_cft_parameters(base_cft, theta, parameter_names)
    _enzyme_continuous_transition!(
        state, cft, global_parameters, climate, day, observable, layer_depth,
        irrigation, false, pathway,
    )
    return _daily_transition_observable(state, observable)
end

include("neural_training.jl")
include("management_adaptation.jl")

end
