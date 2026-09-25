"""Fixed-shape observation contract and optional Enzyme entry points."""

const _AD_TARGETS = (:gpp, :reco, :et)
const _NEURAL_TARGETS = (:gpp, :reco, :et)

"""
    ADSeasonContext(growth_mask, valid_masks, observations, scales; calendar_days=365)

Immutable station-season observation context for the differentiable adapter.
The modeled growth-season mask is authoritative; observations outside that mask
are never counted in the loss. The contract uses the model's no-leap calendar.
"""
struct ADSeasonContext{M, V, O, S, C}
    growth_mask::M
    valid_masks::V
    observations::O
    scales::S
    counts::C
    target_count::Int
    calendar_days::Int
end

"""
    NeuralSeasonContext(growth_mask, valid_masks, observations, scales;
                        weights=nothing, calendar_days=365)

Observation contract for end-to-end NeuralCrop training. GPP, RECO, and ET are
normalized independently; internal soil-water states remain process variables
and are not observation targets.
"""
struct NeuralSeasonContext{M, V, O, S, C, W, T}
    growth_mask::M
    valid_masks::V
    observations::O
    scales::S
    counts::C
    weights::W
    weight_sum::T
    target_count::Int
    calendar_days::Int
end

Adapt.@adapt_structure NeuralSeasonContext

function NeuralSeasonContext(
    growth_mask::AbstractVector{Bool},
    valid_masks::NamedTuple,
    observations::NamedTuple,
    scales::NamedTuple;
    weights = nothing,
    calendar_days::Integer = 365,
)
    calendar_days == 365 || throw(ArgumentError(
        "NeuralSeasonContext requires the 365-day no-leap calendar",
    ))
    days = length(growth_mask)
    days > 0 || throw(ArgumentError("growth_mask must not be empty"))
    valid_values = ntuple(length(_NEURAL_TARGETS)) do index
        name = _NEURAL_TARGETS[index]
        hasproperty(valid_masks, name) || throw(ArgumentError(
            "valid_masks is missing target $name",
        ))
        mask = BitVector(getproperty(valid_masks, name))
        length(mask) == days || throw(DimensionMismatch(
            "valid mask $name has $(length(mask)) rows, expected $days",
        ))
        mask
    end
    observation_values = ntuple(length(_NEURAL_TARGETS)) do index
        name = _NEURAL_TARGETS[index]
        hasproperty(observations, name) || throw(ArgumentError(
            "observations is missing target $name",
        ))
        values = getproperty(observations, name)
        values isa AbstractVector || throw(ArgumentError(
            "observations.$name must be a vector",
        ))
        length(values) == days || throw(DimensionMismatch(
            "observations.$name has $(length(values)) rows, expected $days",
        ))
        values
    end
    scale_values = ntuple(length(_NEURAL_TARGETS)) do index
        name = _NEURAL_TARGETS[index]
        hasproperty(scales, name) || throw(ArgumentError("scales is missing target $name"))
        scale = getproperty(scales, name)
        scale isa Real && isfinite(scale) && scale > 0 || throw(ArgumentError(
            "scales.$name must be finite and positive",
        ))
        scale
    end
    weight_values = if weights === nothing
        ntuple(_ -> 1.0, length(_NEURAL_TARGETS))
    else
        ntuple(length(_NEURAL_TARGETS)) do index
            name = _NEURAL_TARGETS[index]
            hasproperty(weights, name) || throw(ArgumentError(
                "weights is missing target $name",
            ))
            weight = getproperty(weights, name)
            weight isa Real && isfinite(weight) && weight >= 0 || throw(ArgumentError(
                "weights.$name must be finite and non-negative",
            ))
            weight
        end
    end

    growth = BitVector(growth_mask)
    valid_named = NamedTuple{_NEURAL_TARGETS}(valid_values)
    observation_named = NamedTuple{_NEURAL_TARGETS}(observation_values)
    for target_index in eachindex(observation_values)
        name = _NEURAL_TARGETS[target_index]
        valid = valid_values[target_index]
        values = observation_values[target_index]
        for day in eachindex(growth)
            if growth[day] && valid[day] && !isfinite(values[day])
                throw(ArgumentError(
                    "observations.$name contains a non-finite value at a valid day",
                ))
            end
        end
    end
    scale_named = NamedTuple{_NEURAL_TARGETS}(scale_values)
    weight_named = NamedTuple{_NEURAL_TARGETS}(weight_values)
    counts = NamedTuple{_NEURAL_TARGETS}(
        ntuple(
            index -> count(identity, growth .& getfield(valid_named, _NEURAL_TARGETS[index])),
            length(_NEURAL_TARGETS),
        ),
    )
    active = ntuple(
        index -> counts[index] > 0 && weight_values[index] > 0,
        length(_NEURAL_TARGETS),
    )
    target_count = count(identity, active)
    target_count > 0 || throw(ArgumentError(
        "at least one positively weighted target must have valid observations",
    ))
    weight_sum = sum(
        weight_values[index] for index in eachindex(weight_values) if active[index]
    )
    return NeuralSeasonContext(
        growth,
        valid_named,
        observation_named,
        scale_named,
        counts,
        weight_named,
        weight_sum,
        target_count,
        Int(calendar_days),
    )
end

"""
    ManagementAdaptationContext(fertilizer_days, fertilizer_fractions; nitrogen_cost=0)

Fixed discrete management schedule for the optional Enzyme management-adaptation
path. `fertilizer_days` are preselected calendar days and
`fertilizer_fractions` must be non-negative fractions summing to one. The
continuous optimization variable is the total mineral nitrogen input; event
days themselves remain fixed outer-loop choices.
"""
struct ManagementAdaptationContext{D, F, T}
    fertilizer_days::D
    fertilizer_fractions::F
    nitrogen_cost::T
end

function ManagementAdaptationContext(
    fertilizer_days::Tuple,
    fertilizer_fractions::Tuple;
    nitrogen_cost::Real = 0,
)
    length(fertilizer_days) == length(fertilizer_fractions) || throw(DimensionMismatch(
        "fertilizer_days and fertilizer_fractions must have identical lengths",
    ))
    isempty(fertilizer_days) && throw(ArgumentError("at least one fertilizer event is required"))
    all(day -> day isa Integer && day > 0, fertilizer_days) || throw(ArgumentError(
        "fertilizer_days must contain positive integers",
    ))
    issorted(fertilizer_days) && allunique(fertilizer_days) || throw(ArgumentError(
        "fertilizer_days must be strictly increasing",
    ))
    all(fraction -> fraction isa Real && isfinite(fraction) && fraction >= 0, fertilizer_fractions) ||
        throw(ArgumentError("fertilizer_fractions must be finite and non-negative"))
    sum(fertilizer_fractions) ≈ one(sum(fertilizer_fractions)) || throw(ArgumentError(
        "fertilizer_fractions must sum to one",
    ))
    nitrogen_cost isa Real && isfinite(nitrogen_cost) && nitrogen_cost >= 0 || throw(ArgumentError(
        "nitrogen_cost must be finite and non-negative",
    ))
    return ManagementAdaptationContext(fertilizer_days, fertilizer_fractions, nitrogen_cost)
end

function ADSeasonContext(
    growth_mask::AbstractVector{Bool},
    valid_masks::NamedTuple,
    observations::NamedTuple,
    scales::NamedTuple;
    calendar_days::Integer = 365,
)
    calendar_days == 365 || throw(ArgumentError(
        "ADSeasonContext requires the 365-day no-leap calendar",
    ))
    days = length(growth_mask)
    days > 0 || throw(ArgumentError("growth_mask must not be empty"))

    valid_values = ntuple(length(_AD_TARGETS)) do index
        name = _AD_TARGETS[index]
        hasproperty(valid_masks, name) || throw(ArgumentError(
            "valid_masks is missing target $name",
        ))
        mask = BitVector(getproperty(valid_masks, name))
        length(mask) == days || throw(DimensionMismatch(
            "valid mask $name has $(length(mask)) rows, expected $days",
        ))
        mask
    end
    observation_values = ntuple(length(_AD_TARGETS)) do index
        name = _AD_TARGETS[index]
        hasproperty(observations, name) || throw(ArgumentError(
            "observations is missing target $name",
        ))
        values = getproperty(observations, name)
        values isa AbstractVector || throw(ArgumentError(
            "observations.$name must be a vector",
        ))
        length(values) == days || throw(DimensionMismatch(
            "observations.$name has $(length(values)) rows, expected $days",
        ))
        values
    end
    scale_values = ntuple(length(_AD_TARGETS)) do index
        name = _AD_TARGETS[index]
        hasproperty(scales, name) || throw(ArgumentError(
            "scales is missing target $name",
        ))
        scale = getproperty(scales, name)
        scale isa Real && isfinite(scale) && scale > 0 || throw(ArgumentError(
            "scales.$name must be finite and positive",
        ))
        scale
    end

    growth = BitVector(growth_mask)
    for index in eachindex(observation_values)
        valid = valid_values[index]
        observations_for_target = observation_values[index]
        for day in eachindex(growth)
            if growth[day] && valid[day] && !isfinite(observations_for_target[day])
                throw(ArgumentError(
                    "observations.$(_AD_TARGETS[index]) contains a non-finite " *
                    "value at a valid growth-season day",
                ))
            end
        end
    end

    valid_named = NamedTuple{_AD_TARGETS}(valid_values)
    observation_named = NamedTuple{_AD_TARGETS}(observation_values)
    scale_named = NamedTuple{_AD_TARGETS}(scale_values)
    counts = NamedTuple{_AD_TARGETS}(
        ntuple(
            index -> count(identity, growth .& getfield(valid_named, _AD_TARGETS[index])),
            length(_AD_TARGETS),
        ),
    )
    target_count = count(value -> value > 0, values(counts))
    target_count > 0 || throw(ArgumentError(
        "at least one target must have valid observations in the growth season",
    ))
    return ADSeasonContext(
        growth,
        valid_named,
        observation_named,
        scale_named,
        counts,
        target_count,
        Int(calendar_days),
    )
end

"""Fallback raised when Enzyme has not been loaded as an optional extension."""
function enzyme_forward_directional(args...; kwargs...)
    throw(ArgumentError(
        "Enzyme is not loaded; add Enzyme to the active environment before " *
        "calling enzyme_forward_directional",
    ))
end

"""Fallback raised when Enzyme has not been loaded as an optional extension."""
function enzyme_forward_gradient(args...; kwargs...)
    throw(ArgumentError(
        "Enzyme is not loaded; add Enzyme to the active environment before " *
        "calling enzyme_forward_gradient",
    ))
end

"""Fallback raised when Enzyme has not been loaded as an optional extension."""
function enzyme_zero_tangent(args...; kwargs...)
    throw(ArgumentError(
        "Enzyme is not loaded; add Enzyme to the active environment before " *
        "calling enzyme_zero_tangent",
    ))
end

"""Fallback raised when Enzyme has not been loaded as an optional extension."""
function enzyme_prepare_daily_state!(args...; kwargs...)
    throw(ArgumentError(
        "Enzyme is not loaded; add Enzyme to the active environment before " *
        "calling enzyme_prepare_daily_state!",
    ))
end

"""Fallback raised when Enzyme has not been loaded as an optional extension."""
function enzyme_daily_transition_objective(args...; kwargs...)
    throw(ArgumentError(
        "Enzyme is not loaded; add Enzyme to the active environment before " *
        "calling enzyme_daily_transition_objective",
    ))
end

"""Fallback raised when Enzyme has not been loaded as an optional extension."""
function enzyme_seasonal_loss(args...; kwargs...)
    throw(ArgumentError(
        "Enzyme is not loaded; add Enzyme to the active environment before " *
        "calling enzyme_seasonal_loss",
    ))
end

"""Fallback raised when Enzyme is not loaded as an optional extension."""
function enzyme_seasonal_gradient_blockwise(args...; kwargs...)
    throw(ArgumentError(
        "Enzyme is not loaded; add Enzyme to the active environment before " *
        "calling enzyme_seasonal_gradient_blockwise",
    ))
end

"""Fallback raised when Enzyme has not been loaded as an optional extension."""
function enzyme_seasonal_soil_loss(args...; kwargs...)
    throw(ArgumentError(
        "Enzyme is not loaded; add Enzyme to the active environment before " *
        "calling enzyme_seasonal_soil_loss",
    ))
end

"""Fallback raised when Enzyme has not been loaded as an optional extension."""
function enzyme_seasonal_soil_gradient_blockwise(args...; kwargs...)
    throw(ArgumentError(
        "Enzyme is not loaded; add Enzyme to the active environment before " *
        "calling enzyme_seasonal_soil_gradient_blockwise",
    ))
end

"""Fallback raised when Enzyme has not been loaded as an optional extension."""
function enzyme_seasonal_joint_loss(args...; kwargs...)
    throw(ArgumentError(
        "Enzyme is not loaded; add Enzyme to the active environment before " *
        "calling enzyme_seasonal_joint_loss",
    ))
end

"""Fallback raised when Enzyme has not been loaded as an optional extension."""
function enzyme_seasonal_joint_gradient_blockwise(args...; kwargs...)
    throw(ArgumentError(
        "Enzyme is not loaded; add Enzyme to the active environment before " *
        "calling enzyme_seasonal_joint_gradient_blockwise",
    ))
end

"""Fallback raised when Enzyme has not been loaded as an optional extension."""
function enzyme_neural_seasonal_loss(args...; kwargs...)
    throw(ArgumentError(
        "Enzyme is not loaded; add Enzyme to the active environment before " *
        "calling enzyme_neural_seasonal_loss",
    ))
end

function enzyme_neural_seasonal_loss_buffer!(args...; kwargs...)
    throw(ArgumentError(
        "load Enzyme.jl before calling enzyme_neural_seasonal_loss_buffer!",
    ))
end

"""Fallback raised when Enzyme has not been loaded as an optional extension."""
function enzyme_neural_seasonal_gradient_blockwise(args...; kwargs...)
    throw(ArgumentError(
        "Enzyme is not loaded; add Enzyme to the active environment before " *
        "calling enzyme_neural_seasonal_gradient_blockwise",
    ))
end

"""Fallback raised when Enzyme has not been loaded as an optional extension."""
function enzyme_neural_predictions(args...; kwargs...)
    throw(ArgumentError(
        "Enzyme is not loaded; add Enzyme to the active environment before " *
        "calling enzyme_neural_predictions",
    ))
end

"""Fallback raised when Enzyme has not been loaded as an optional extension."""
function enzyme_management_yield_loss(args...; kwargs...)
    throw(ArgumentError(
        "Enzyme is not loaded; add Enzyme to the active environment before " *
        "calling enzyme_management_yield_loss",
    ))
end

"""Fallback raised when Enzyme has not been loaded as an optional extension."""
function enzyme_management_yield_split_loss(args...; kwargs...)
    throw(ArgumentError(
        "Enzyme is not loaded; add Enzyme to the active environment before " *
        "calling enzyme_management_yield_split_loss",
    ))
end

"""Fallback raised when Enzyme has not been loaded as an optional extension."""
function enzyme_joint_adaptation_yield_loss(args...; kwargs...)
    throw(ArgumentError(
        "Enzyme is not loaded; add Enzyme to the active environment before " *
        "calling enzyme_joint_adaptation_yield_loss",
    ))
end
