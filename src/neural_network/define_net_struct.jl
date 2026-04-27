struct MLP{M <: Lux.AbstractLuxLayer, K} <:Lux.Lux.AbstractLuxWrapperLayer{:model}
    model::M 
    kwargs::K
end

function MLP(model::Lux.AbstractLuxLayer;  kwargs...)
    return MLP{typeof(model), typeof(kwargs)}(model, kwargs)
end

struct NODE{M <: Lux.AbstractLuxLayer, So, T, K} <:Lux.Lux.AbstractLuxWrapperLayer{:model}
    model::M
    solver::So
    tspan::T
    kwargs::K
end

function NODE(model::Lux.AbstractLuxLayer; solver=nothing, tspan=(0.0f0, 1.0f0), kwargs...)
    return NODE{typeof(model), typeof(solver), typeof(tspan), typeof(kwargs)}(model, solver, tspan, kwargs)
end