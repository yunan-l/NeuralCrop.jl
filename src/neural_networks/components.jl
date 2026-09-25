"""A two-hidden-layer MLP view into one shared, flat parameter vector."""
struct MLPLayout{NInput, NHidden, NOutput}
    first_parameter::Int
end

@inline function _mlp_parameter_count(::MLPLayout{NInput, NHidden, NOutput}) where {
    NInput, NHidden, NOutput,
}
    return NHidden * NInput + NHidden + NHidden * NHidden + NHidden +
           NOutput * NHidden + NOutput
end

@inline _next_parameter(layout::MLPLayout) =
    layout.first_parameter + _mlp_parameter_count(layout)

"""
    NeuralComponents(; kwargs...)

Select the retained NeuralCrop components. A disabled component falls back to its
ordinary NeuralCrop process implementation, which is used for controlled
ablation runs. GPP can be predicted directly or used as a bounded
multiplicative correction to process GPP. The older lambda and Vcmax networks
are retained only so archived checkpoints remain describable. Carbon
allocation, soil decomposition, snowmelt, and soil evaporation are likewise
retained but disabled in the current experiments.
"""
Base.@kwdef struct NeuralComponents
    gpp::Bool = true
    gpp_residual::Bool = false
    lambda::Bool = false
    vcmax::Bool = false
    allocation::Bool = false
    decomposition::Bool = false
    respiration::Bool = true
    snowmelt::Bool = false
    transpiration::Bool = true
    evaporation::Bool = false
end

"""
    NeuralCropLayout(; hidden=64, components=NeuralComponents())

Flat parameter layout for the retained NeuralCrop networks. The production
contract intentionally fixes two hidden layers with 64 units each. Inputs are
limited to the variables described in the manuscript:

* GPP: day length, correctly scaled APAR, fPAR, LAI, leaf nitrogen, air
  temperature, atmospheric CO2, and mean top-three-layer soil moisture;
* lambda and Vcmax: legacy networks retained only for archived checkpoints;
* allocation: NPP, LAI, leaf nitrogen, soil moisture, and four carbon pools;
* decomposition: temperature and moisture;
* respiration: GPP, leaf respiration, temperature, and crop carbon pools;
* snowmelt: air temperature, snowpack, and precipitation;
* transpiration: day length, equilibrium PET, canopy state, root carbon, and
  root-weighted soil moisture; and
* evaporation: equilibrium PET, vegetation/wetness/cover, and top-three water.

The lambda, Vcmax, allocation, decomposition, snowmelt, and soil-evaporation
blocks are inactive by default. Their native NeuralCrop process equations are
used in the formal experiments.
"""
struct NeuralCropLayout{G, L, V, A, D, R, S, T, E, F}
    gpp::G
    lambda::L
    vcmax::V
    allocation::A
    decomposition::D
    respiration::R
    snowmelt::S
    transpiration::T
    evaporation::E
    components::F
end

function NeuralCropLayout(;
    hidden::Integer = 64,
    components::NeuralComponents = NeuralComponents(),
)
    hidden == 64 || throw(ArgumentError(
        "NeuralCrop's frozen architecture requires 64 hidden units",
    ))
    components.gpp_residual && !components.gpp && throw(ArgumentError(
        "gpp_residual requires the GPP network to be enabled",
    ))
    first = 1
    gpp = MLPLayout{8, 64, 1}(first)
    first = _next_parameter(gpp)
    lambda = MLPLayout{3, 64, 1}(first)
    first = _next_parameter(lambda)
    vcmax = MLPLayout{4, 64, 1}(first)
    first = _next_parameter(vcmax)
    allocation = MLPLayout{8, 64, 4}(first)
    first = _next_parameter(allocation)
    decomposition = MLPLayout{2, 64, 1}(first)
    first = _next_parameter(decomposition)
    respiration = MLPLayout{7, 64, 1}(first)
    first = _next_parameter(respiration)
    snowmelt = MLPLayout{3, 64, 1}(first)
    first = _next_parameter(snowmelt)
    transpiration = MLPLayout{7, 64, 1}(first)
    first = _next_parameter(transpiration)
    evaporation = MLPLayout{7, 64, 4}(first)
    return NeuralCropLayout(
        gpp, lambda, vcmax, allocation, decomposition, respiration, snowmelt,
        transpiration, evaporation, components,
    )
end

"""Number of stored scalar parameters, including inactive retained blocks."""
@inline neural_parameter_count(layout::NeuralCropLayout) =
    _next_parameter(layout.evaporation) - 1

"""Named ranges for checkpointing and component-wise ablation diagnostics."""
function neural_parameter_ranges(layout::NeuralCropLayout)
    component_range(component) =
        component.first_parameter:(_next_parameter(component) - 1)
    return (
        gpp = component_range(layout.gpp),
        lambda = component_range(layout.lambda),
        vcmax = component_range(layout.vcmax),
        allocation = component_range(layout.allocation),
        decomposition = component_range(layout.decomposition),
        respiration = component_range(layout.respiration),
        snowmelt = component_range(layout.snowmelt),
        transpiration = component_range(layout.transpiration),
        evaporation = component_range(layout.evaporation),
    )
end

"""Number of parameters enabled by the component selection."""
function neural_trainable_parameter_count(layout::NeuralCropLayout)
    ranges = neural_parameter_ranges(layout)
    components = layout.components
    return (components.gpp ? length(ranges.gpp) : 0) +
           (components.lambda ? length(ranges.lambda) : 0) +
           (components.vcmax ? length(ranges.vcmax) : 0) +
           (components.allocation ? length(ranges.allocation) : 0) +
           (components.decomposition ? length(ranges.decomposition) : 0) +
           (components.respiration ? length(ranges.respiration) : 0) +
           (components.snowmelt ? length(ranges.snowmelt) : 0) +
           (components.transpiration ? length(ranges.transpiration) : 0) +
           (components.evaporation ? length(ranges.evaporation) : 0)
end

@inline _signed_scale(value::T, scale::T) where {T <: AbstractFloat} =
    value / (abs(value) + scale)

@inline _positive_scale(value::T, scale::T) where {T <: AbstractFloat} =
    max(value, zero(T)) / (max(value, zero(T)) + scale)

@inline _logistic(value::T) where {T <: AbstractFloat} =
    value >= zero(T) ? inv(one(T) + exp(-value)) : exp(value) / (one(T) + exp(value))

@inline function _mlp_forward(
    theta::AbstractVector{T},
    layout::MLPLayout{NInput, NHidden, NOutput},
    inputs::NTuple{NInput, T},
) where {T <: AbstractFloat, NInput, NHidden, NOutput}
    first = layout.first_parameter
    w1 = first
    b1 = w1 + NHidden * NInput
    w2 = b1 + NHidden
    b2 = w2 + NHidden * NHidden
    w3 = b2 + NHidden
    b3 = w3 + NOutput * NHidden

    hidden1 = ntuple(Val(NHidden)) do row
        value = theta[b1 + row - 1]
        @inbounds for column in 1:NInput
            value += theta[w1 + (row - 1) * NInput + column - 1] * inputs[column]
        end
        tanh(value)
    end
    hidden2 = ntuple(Val(NHidden)) do row
        value = theta[b2 + row - 1]
        @inbounds for column in 1:NHidden
            value += theta[w2 + (row - 1) * NHidden + column - 1] * hidden1[column]
        end
        tanh(value)
    end
    return ntuple(Val(NOutput)) do row
        value = theta[b3 + row - 1]
        @inbounds for column in 1:NHidden
            value += theta[w3 + (row - 1) * NHidden + column - 1] * hidden2[column]
        end
        value
    end
end

function _initialize_mlp!(
    rng::AbstractRNG,
    theta::AbstractVector{T},
    layout::MLPLayout{NInput, NHidden, NOutput},
    output_biases::NTuple{NOutput, T},
    ; zero_output_weights::Bool = false,
) where {T <: AbstractFloat, NInput, NHidden, NOutput}
    first = layout.first_parameter
    w1 = first
    b1 = w1 + NHidden * NInput
    w2 = b1 + NHidden
    b2 = w2 + NHidden * NHidden
    w3 = b2 + NHidden
    b3 = w3 + NOutput * NHidden
    first_scale = sqrt(T(2) / T(NInput + NHidden))
    hidden_scale = sqrt(T(2) / T(2 * NHidden))
    output_scale = T(1e-3)
    @inbounds for index in w1:(b1 - 1)
        theta[index] = first_scale * randn(rng, T)
    end
    fill!(view(theta, b1:(w2 - 1)), zero(T))
    @inbounds for index in w2:(b2 - 1)
        theta[index] = hidden_scale * randn(rng, T)
    end
    fill!(view(theta, b2:(w3 - 1)), zero(T))
    @inbounds for index in w3:(b3 - 1)
        theta[index] = zero_output_weights ? zero(T) : output_scale * randn(rng, T)
    end
    @inbounds for output in 1:NOutput
        theta[b3 + output - 1] = output_biases[output]
    end
    return theta
end

"""
    initialize_neural_parameters(layout; T=Float64, seed=20260921)

Create deterministic Xavier-style parameters. Output biases start every
component in a physically moderate regime without pretraining against process
model output.
"""
function initialize_neural_parameters(
    layout::NeuralCropLayout;
    T::Type{<:AbstractFloat} = Float64,
    seed::Integer = 20260921,
)
    rng = Xoshiro(seed)
    theta = zeros(T, neural_parameter_count(layout))
    if layout.components.gpp_residual
        _initialize_mlp!(
            rng, theta, layout.gpp, (zero(T),); zero_output_weights = true,
        )
    else
        _initialize_mlp!(rng, theta, layout.gpp, (log(T(0.2) / T(0.8)),))
    end
    _initialize_mlp!(rng, theta, layout.lambda, (T(1.5),))
    _initialize_mlp!(rng, theta, layout.vcmax, (zero(T),))
    _initialize_mlp!(rng, theta, layout.allocation,
        (log(T(0.3)), log(T(0.3)), log(T(0.2)), log(T(0.2))))
    _initialize_mlp!(rng, theta, layout.decomposition, (T(-1.4),))
    _initialize_mlp!(
        rng, theta, layout.respiration, (zero(T),); zero_output_weights = true,
    )
    _initialize_mlp!(rng, theta, layout.snowmelt, (T(-2.2),))
    _initialize_mlp!(
        rng, theta, layout.transpiration, (zero(T),); zero_output_weights = true,
    )
    _initialize_mlp!(rng, theta, layout.evaporation,
        (T(-0.85), zero(T), zero(T), zero(T)))
    return theta
end

"""Bounded multiplicative correction to process GPP in the interval `[0, 2]`."""
@inline function neural_gpp_multiplier(
    theta::AbstractVector{T}, layout::NeuralCropLayout,
    daylength::T, apar::T, fpar::T, lai::T, leaf_nitrogen::T,
    air_temperature::T, co2::T, top3_moisture::T,
) where {T <: AbstractFloat}
    inputs = (
        daylength / T(12) - one(T),
        T(2) * _positive_scale(apar, T(1e7)) - one(T),
        T(2) * clamp(fpar, zero(T), one(T)) - one(T),
        T(2) * _positive_scale(lai, T(5)) - one(T),
        T(2) * _positive_scale(leaf_nitrogen, T(5)) - one(T),
        _signed_scale(air_temperature, T(25)),
        clamp((co2 - T(400)) / T(200), -one(T), one(T)),
        T(2) * clamp(top3_moisture, zero(T), one(T)) - one(T),
    )
    return one(T) + tanh(_mlp_forward(theta, layout.gpp, inputs)[1])
end

"""Direct neural prediction of daily crop GPP in g C m^-2 day^-1."""
@inline function neural_gpp(
    theta::AbstractVector{T}, layout::NeuralCropLayout,
    daylength::T, apar::T, fpar::T, lai::T, leaf_nitrogen::T,
    air_temperature::T, co2::T, top3_moisture::T,
) where {T <: AbstractFloat}
    inputs = (
        daylength / T(12) - one(T),
        T(2) * _positive_scale(apar, T(1e7)) - one(T),
        T(2) * clamp(fpar, zero(T), one(T)) - one(T),
        T(2) * _positive_scale(lai, T(5)) - one(T),
        T(2) * _positive_scale(leaf_nitrogen, T(5)) - one(T),
        _signed_scale(air_temperature, T(25)),
        clamp((co2 - T(400)) / T(200), -one(T), one(T)),
        T(2) * clamp(top3_moisture, zero(T), one(T)) - one(T),
    )
    return T(50) * _logistic(_mlp_forward(theta, layout.gpp, inputs)[1])
end

"""Neural replacement for the daily intercellular-to-ambient CO2 ratio."""
@inline function neural_lambda(
    theta::AbstractVector{T}, layout::NeuralCropLayout,
    daylength::T, air_temperature::T, top3_moisture::T,
) where {T <: AbstractFloat}
    inputs = (
        daylength / T(12) - one(T),
        _signed_scale(air_temperature, T(25)),
        T(2) * clamp(top3_moisture, zero(T), one(T)) - one(T),
    )
    raw = _mlp_forward(theta, layout.lambda, inputs)[1]
    return _logistic(raw)
end

"""Neural replacement for daily maximum Rubisco capacity."""
@inline function neural_vcmax(
    theta::AbstractVector{T}, layout::NeuralCropLayout,
    daylength::T, apar::T, leaf_nitrogen::T, air_temperature::T,
) where {T <: AbstractFloat}
    inputs = (
        daylength / T(12) - one(T),
        _positive_scale(apar, T(20)),
        _positive_scale(leaf_nitrogen, T(5)),
        _signed_scale(air_temperature, T(25)),
    )
    raw = _mlp_forward(theta, layout.vcmax, inputs)[1]
    return T(200) * _logistic(raw)
end

"""Mass-conserving leaf/root/storage/pool fractions for the retained inactive block."""
@inline function neural_allocation_fractions(
    theta::AbstractVector{T}, layout::NeuralCropLayout,
    npp::T, lai::T, leaf_nitrogen::T, top3_moisture::T,
    leaf_carbon::T, root_carbon::T, storage_carbon::T, pool_carbon::T,
) where {T <: AbstractFloat}
    inputs = (
        _signed_scale(npp, T(10)),
        _positive_scale(lai, T(5)),
        _positive_scale(leaf_nitrogen, T(5)),
        T(2) * clamp(top3_moisture, zero(T), one(T)) - one(T),
        _positive_scale(leaf_carbon, T(500)),
        _positive_scale(root_carbon, T(500)),
        _positive_scale(storage_carbon, T(500)),
        _positive_scale(pool_carbon, T(500)),
    )
    logits = _mlp_forward(theta, layout.allocation, inputs)
    maximum = max(max(logits[1], logits[2]), max(logits[3], logits[4]))
    weights = ntuple(Val(4)) do index
        exp(logits[index] - maximum)
    end
    total = weights[1] + weights[2] + weights[3] + weights[4]
    return ntuple(Val(4)) do index
        weights[index] / total
    end
end

"""Bounded neural environmental response shared by carbon and nitrogen turnover."""
@inline function neural_decomposition_response(
    theta::AbstractVector{T}, layout::NeuralCropLayout,
    temperature::T, moisture::T,
) where {T <: AbstractFloat}
    temperature >= T(-15) || return zero(T)
    inputs = (
        _signed_scale(temperature, T(25)),
        T(2) * clamp(moisture, zero(T), one(T)) - one(T),
    )
    return _logistic(_mlp_forward(theta, layout.decomposition, inputs)[1])
end

"""Bounded multiplier around process crop respiration; initialization is one."""
@inline function neural_crop_respiration_multiplier(
    theta::AbstractVector{T}, layout::NeuralCropLayout,
    gross_assimilation::T, leaf_respiration::T,
    air_temperature::T, soil_temperature::T,
    root_carbon::T, storage_carbon::T, pool_carbon::T,
) where {T <: AbstractFloat}
    inputs = (
        _positive_scale(gross_assimilation, T(20)),
        _positive_scale(leaf_respiration, T(5)),
        _signed_scale(air_temperature, T(25)),
        _signed_scale(soil_temperature, T(25)),
        _positive_scale(root_carbon, T(500)),
        _positive_scale(storage_carbon, T(500)),
        _positive_scale(pool_carbon, T(500)),
    )
    return T(2) * _logistic(_mlp_forward(theta, layout.respiration, inputs)[1])
end

"""Bounded snowmelt flux retained as an optional, default-off component."""
@inline function neural_snow_melt(
    theta::AbstractVector{T}, layout::NeuralCropLayout,
    air_temperature::T, snowpack::T, precipitation::T,
) where {T <: AbstractFloat}
    (air_temperature > zero(T) && snowpack > zero(T)) || return zero(T)
    inputs = (
        _signed_scale(air_temperature, T(10)),
        _positive_scale(snowpack, T(100)),
        _positive_scale(precipitation, T(20)),
    )
    fraction = _logistic(_mlp_forward(theta, layout.snowmelt, inputs)[1])
    return min(snowpack, snowpack * fraction)
end

"""Bounded multiplier around process transpiration; initialization is one."""
@inline function neural_transpiration_multiplier(
    theta::AbstractVector{T}, layout::NeuralCropLayout,
    daylength::T, equilibrium_pet::T, fpar::T, canopy_wet::T,
    canopy_conductance::T, root_carbon::T, root_water::T,
) where {T <: AbstractFloat}
    inputs = (
        daylength / T(12) - one(T),
        _positive_scale(equilibrium_pet, T(5)),
        T(2) * clamp(fpar, zero(T), one(T)) - one(T),
        T(2) * clamp(canopy_wet, zero(T), one(T)) - one(T),
        _positive_scale(canopy_conductance, T(1)),
        _positive_scale(root_carbon, T(500)),
        T(2) * clamp(root_water, zero(T), one(T)) - one(T),
    )
    return T(2) * _logistic(_mlp_forward(theta, layout.transpiration, inputs)[1])
end

"""
    neural_soil_evaporation(...)

Return a bounded fraction of physically available evaporation followed by
three non-negative allocation weights. Only the first three soil layers are
learned; deeper-layer evaporation remains process based.
"""
@inline function neural_soil_evaporation(
    theta::AbstractVector{T}, layout::NeuralCropLayout,
    equilibrium_pet::T, fpar::T, canopy_wet::T, litter_cover::T,
    water1::T, water2::T, water3::T,
) where {T <: AbstractFloat}
    inputs = (
        _positive_scale(equilibrium_pet, T(5)),
        T(2) * clamp(fpar, zero(T), one(T)) - one(T),
        T(2) * clamp(canopy_wet, zero(T), one(T)) - one(T),
        T(2) * clamp(litter_cover, zero(T), one(T)) - one(T),
        _positive_scale(water1, T(40)),
        _positive_scale(water2, T(60)),
        _positive_scale(water3, T(100)),
    )
    outputs = _mlp_forward(theta, layout.evaporation, inputs)
    total_fraction = _logistic(outputs[1])
    maximum = max(outputs[2], max(outputs[3], outputs[4]))
    weights = (
        exp(outputs[2] - maximum),
        exp(outputs[3] - maximum),
        exp(outputs[4] - maximum),
    )
    total = weights[1] + weights[2] + weights[3]
    return (total_fraction, weights[1] / total, weights[2] / total, weights[3] / total)
end
