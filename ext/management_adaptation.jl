"""Fixed-event, continuous management objectives for future-scenario studies."""

@inline function _management_fertilizer_fraction(context, day::Integer, ::Type{T}) where {T}
    fraction = zero(T)
    for index in eachindex(context.fertilizer_days)
        fraction += day == context.fertilizer_days[index] ?
            T(context.fertilizer_fractions[index]) : zero(T)
    end
    return fraction
end

"""
    _apply_management_fertilizer!(state, total_fertilizer, context, day, lpjmlparams)

Apply one fixed-date mineral-nitrogen pulse to the top soil layer. This is the
continuous counterpart of the prescribed-fertilizer split: nitrate and
ammonium fractions remain those in `LPJmLParams`, while pulse dates and pulse
fractions are immutable outer inputs. The method intentionally supports one
cell only; global optimization will batch independent cell objectives later.
"""
@inline function _apply_management_fertilizer!(
    state::NeuralCrop.ModelState,
    total_fertilizer::T,
    context::NeuralCrop.ManagementAdaptationContext,
    day::Integer,
    lpjmlparams,
) where {T <: AbstractFloat}
    nitrogen = NeuralCrop.soil_nitrogen_prognostic(state)
    size(nitrogen.nitrate, 2) == 1 || throw(ArgumentError(
        "enzyme_management_yield_loss currently supports one cell per objective",
    ))
    pulse = total_fertilizer * _management_fertilizer_fraction(context, day, T)
    nitrogen.nitrate[1, 1] += pulse * T(lpjmlparams.nfert_no3_frac)
    nitrogen.ammonium[1, 1] += pulse * (one(T) - T(lpjmlparams.nfert_no3_frac))
    return nothing
end

@inline function _management_split_fraction(
    context::NeuralCrop.ManagementAdaptationContext,
    day::Integer,
    first_fraction::T,
) where {T <: AbstractFloat}
    first_day, second_day = context.fertilizer_days
    return day == first_day ? first_fraction :
        day == second_day ? one(T) - first_fraction : zero(T)
end

@inline function _apply_management_split_fertilizer!(
    state::NeuralCrop.ModelState,
    total_fertilizer::T,
    first_fraction::T,
    context::NeuralCrop.ManagementAdaptationContext,
    day::Integer,
    lpjmlparams,
) where {T <: AbstractFloat}
    nitrogen = NeuralCrop.soil_nitrogen_prognostic(state)
    size(nitrogen.nitrate, 2) == 1 || throw(ArgumentError(
        "enzyme_management_yield_split_loss currently supports one cell per objective",
    ))
    fraction = _management_split_fraction(context, day, first_fraction)
    pulse = total_fertilizer * fraction
    nitrogen.nitrate[1, 1] += pulse * T(lpjmlparams.nfert_no3_frac)
    nitrogen.ammonium[1, 1] += pulse * (one(T) - T(lpjmlparams.nfert_no3_frac))
    return nothing
end

"""
    NeuralCrop.enzyme_management_yield_loss(
        theta, state, cft, global_parameters, climate, days, layer_depth, context;
        irrigation=false,
    )

Return the negative fixed-event management benefit for a single C3 cell. The
sole entry in `theta` is total mineral fertilizer (gN m⁻² crop⁻¹), distributed
over immutable event days from `context`. The objective is negative
pre-harvest storage-organ carbon plus nitrogen cost, so minimizing it maximizes
storage carbon net of fertilizer input. Cultivation, harvest, crop failure,
sowing-date selection, and fertilization-event timing are deliberately outside
this differentiable path.
"""
function NeuralCrop.enzyme_management_yield_loss(
    theta::AbstractVector{T},
    state::NeuralCrop.ModelState,
    cft::NeuralCrop.CFTParameters,
    global_parameters,
    climate,
    days::AbstractUnitRange{<:Integer},
    layer_depth,
    context::NeuralCrop.ManagementAdaptationContext;
    irrigation::Bool = false,
) where {T <: AbstractFloat}
    length(theta) == 1 || throw(DimensionMismatch(
        "management adaptation expects theta = [total_fertilizer]",
    ))
    cft.path == 1 || throw(ArgumentError(
        "enzyme_management_yield_loss currently supports C3 CFTs only",
    ))

    model_parameters = _enzyme_model_parameters(T, global_parameters)
    isempty(days) && throw(ArgumentError("management adaptation requires at least one day"))
    day = first(days)
    final_day = last(days)
    while day <= final_day
        _apply_management_fertilizer!(state, theta[1], context, day, model_parameters.lpjml)
        _enzyme_continuous_transition!(
            state, cft, model_parameters, climate, day, :gpp, layer_depth,
            irrigation, false, Val(:C3),
        )
        day += 1
    end
    storage_carbon = NeuralCrop.crop_prognostic(state).carbon.storage[1]
    return -(storage_carbon - T(context.nitrogen_cost) * theta[1])
end

"""
    NeuralCrop.enzyme_management_yield_split_loss(
        theta, state, cft, global_parameters, climate, days, layer_depth, context;
        irrigation=false,
    )

Return the negative fixed-event management benefit for a single C3 cell with
`theta = [total_fertilizer, first_event_fraction]`. The context must contain
exactly two immutable fertilizer dates. The first pulse is `F * r` and the
second is `F * (1-r)`; callers are responsible for constraining `r` to `[0, 1]`.
The fixed fractions stored in `context` remain the contract for the original
one-variable objective and are not used by this two-variable objective.
"""
function NeuralCrop.enzyme_management_yield_split_loss(
    theta::AbstractVector{T},
    state::NeuralCrop.ModelState,
    cft::NeuralCrop.CFTParameters,
    global_parameters,
    climate,
    days::AbstractUnitRange{<:Integer},
    layer_depth,
    context::NeuralCrop.ManagementAdaptationContext;
    irrigation::Bool = false,
) where {T <: AbstractFloat}
    length(theta) == 2 || throw(DimensionMismatch(
        "split management adaptation expects theta = [total_fertilizer, first_event_fraction]",
    ))
    length(context.fertilizer_days) == 2 || throw(ArgumentError(
        "split management adaptation requires exactly two fertilizer event days",
    ))
    cft.path == 1 || throw(ArgumentError(
        "enzyme_management_yield_split_loss currently supports C3 CFTs only",
    ))

    model_parameters = _enzyme_model_parameters(T, global_parameters)
    isempty(days) && throw(ArgumentError("management adaptation requires at least one day"))
    day = first(days)
    final_day = last(days)
    while day <= final_day
        _apply_management_split_fertilizer!(
            state, theta[1], theta[2], context, day, model_parameters.lpjml,
        )
        _enzyme_continuous_transition!(
            state, cft, model_parameters, climate, day, :gpp, layer_depth,
            irrigation, false, Val(:C3),
        )
        day += 1
    end
    storage_carbon = NeuralCrop.crop_prognostic(state).carbon.storage[1]
    return -(storage_carbon - T(context.nitrogen_cost) * theta[1])
end

"""
    NeuralCrop.enzyme_joint_adaptation_yield_loss(
        theta, state, cft, global_parameters, climate, days, layer_depth, context;
        irrigation=false,
    )

Return the fixed-event management-and-cultivar surrogate loss for one C3 cell.
The continuous vector is ordered as
`[total_fertilizer, first_event_fraction, phu_scale,
heat_tolerance_shift, photoperiod_sensitivity_scale,
vernalization_requirement_scale]`.

Only the continuous process trajectory belongs to this objective. Sowing date,
crop establishment, integer `pvd_max`, harvest, and failure remain ordinary
production events. The PHU and vernalization decisions therefore scale the
already-established fixed-event state. Every selected solution must be
materialized and re-evaluated through the ordinary production driver.
"""
function NeuralCrop.enzyme_joint_adaptation_yield_loss(
    theta::AbstractVector{T},
    state::NeuralCrop.ModelState,
    cft::NeuralCrop.CFTParameters,
    global_parameters,
    climate,
    days::AbstractUnitRange{<:Integer},
    layer_depth,
    context::NeuralCrop.ManagementAdaptationContext;
    irrigation::Bool = false,
) where {T <: AbstractFloat}
    length(theta) == 6 || throw(DimensionMismatch(
        "joint adaptation expects [fertilizer, split, PHU scale, heat shift, " *
        "photoperiod scale, vernalization scale]",
    ))
    length(context.fertilizer_days) == 2 || throw(ArgumentError(
        "joint adaptation requires exactly two fertilizer event days",
    ))
    cft.path == 1 || throw(ArgumentError(
        "enzyme_joint_adaptation_yield_loss currently supports C3 CFTs only",
    ))

    adapted_values = T[
        T(cft.temp_photos.high) + theta[4],
        T(cft.temp_co2.high) + theta[4],
        T(cft.psens) * theta[5],
    ]
    adapted_cft = _replace_cft_parameters(
        cft,
        adapted_values,
        (:temp_photos_high, :temp_co2_high, :psens),
    )

    phenology = NeuralCrop.crop_phenology_input(state)
    phenology.phu[1] *= theta[3]
    state.prognostic.climate.V_req[1] *= theta[6]

    model_parameters = _enzyme_model_parameters(T, global_parameters)
    isempty(days) && throw(ArgumentError("joint adaptation requires at least one day"))
    day = first(days)
    final_day = last(days)
    while day <= final_day
        _apply_management_split_fertilizer!(
            state, theta[1], theta[2], context, day, model_parameters.lpjml,
        )
        _enzyme_continuous_transition!(
            state,
            adapted_cft,
            model_parameters,
            climate,
            day,
            :gpp,
            layer_depth,
            irrigation,
            false,
            Val(:C3),
        )
        day += 1
    end
    storage_carbon = NeuralCrop.crop_prognostic(state).carbon.storage[1]
    return -(storage_carbon - T(context.nitrogen_cost) * theta[1])
end
