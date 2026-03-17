struct MLP{M <: Lux.AbstractExplicitLayer, K} <:Lux.AbstractExplicitContainerLayer{(:model,)}
    model::M 
    kwargs::K
end

function MLP(model::Lux.AbstractExplicitLayer;  kwargs...)
    return MLP{typeof(model), typeof(kwargs)}(model, kwargs)
end

struct NODE{M <: Lux.AbstractExplicitLayer, So, T, K} <:Lux.AbstractExplicitContainerLayer{(:model,)}
    model::M
    solver::So
    tspan::T
    kwargs::K
end

function NODE(model::Lux.AbstractExplicitLayer; solver=nothing, tspan=(0.0f0, 1.0f0), kwargs...)
    return NODE{typeof(model), typeof(solver), typeof(tspan), typeof(kwargs)}(model, solver, tspan, kwargs)
end