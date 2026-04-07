# lambda
function neural_lambda(n::MLP, ps, st, input)
    
    st_model = Lux.StatefulLuxLayer{true}(n.model, ps, st)
        
    return vec(st_model(input, ps))
    
end

# vmax
function neural_vmax(n::MLP, ps, st, input)
    
    st_model = Lux.StatefulLuxLayer{true}(n.model, ps, st)
        
    return vec(st_model(input, ps))
    
end

# storage carbon
function neural_stoc(n::NODE, u0, ps, st, input; dt = 1.0f0)
    
    st_model = Lux.StatefulLuxLayer{true}(n.model, ps, st)

    rhs(u, p, t) = st_model(vcat(u/200, input), p)
    
    prob = ODEProblem{false}(ODEFunction{false}(rhs), u0, n.tspan, ps)
    
    return vec(solve(prob, n.solver; dt = dt))
end

# vegetation carbon pools
function neural_allocation(n::NODE, u0, ps, st, input; dt = 1.0f0)
    
    st_model = Lux.StatefulLuxLayer{true}(n.model, ps, st)

    rhs(u, p, t) = st_model(vcat(u/200, input), p)
    
    prob = ODEProblem{false}(ODEFunction{false}(rhs), u0, n.tspan, ps)
    
    return solve(prob, n.solver; dt = dt)
end

# litter carbon pools
function hybrid_litc(n::NODE, u0, ps, st, input, response; dt = 1.0f0)
    
    st_model = Lux.StatefulLuxLayer{true}(n.model, ps, st)
    
    rhs(u, p, t) = -(1.0f0 .- exp.(-response .* st_model(input, p)/365)) .* u
  
    prob = ODEProblem{false}(ODEFunction{false}(rhs), u0, n.tspan, ps) 
    
    return solve(prob, n.solver; dt = dt)
end

# litter nitrogen pools
function hybrid_litn(n::NODE, u0, ps, st, input, response; dt = 1.0f0)
    
    st_model = Lux.StatefulLuxLayer{true}(n.model, ps, st)
    
    rhs(u, p, t) = -(1.0f0 .- exp.(-response .* st_model(input, p)/365)) .* u
  
    prob = ODEProblem{false}(ODEFunction{false}(rhs), u0, n.tspan, ps) 
    
    return solve(prob, n.solver; dt = dt)
end

# soil carbon pools
function hybrid_soilc(n::NODE, u0, ps, st, input, response, A_trans, c_input; dt = 1.0f0)

    st_model = Lux.StatefulLuxLayer{true}(n.model, ps, st)
    
    rhs(u, p, t) = -(1.0f0 .- exp.(-response .* st_model(input, p)/365)) .* u
  
    prob = ODEProblem{false}(ODEFunction{false}(rhs), u0, n.tspan, ps) 
    
    return solve(prob, n.solver; A_trans = A_trans, c_input = c_input, dt = dt)
end

# soil nitrogen pools
function hybrid_soiln(n::NODE, u0, ps, st, input, response, A_trans, c_input; dt = 1.0f0)

    st_model = Lux.StatefulLuxLayer{true}(n.model, ps, st)
    
    rhs(u, p, t) = -(1.0f0 .- exp.(-response .* st_model(input, p)/365)) .* u
  
    prob = ODEProblem{false}(ODEFunction{false}(rhs), u0, n.tspan, ps) 
    
    return solve(prob, n.solver; A_trans = A_trans, c_input = c_input, dt = dt)
end


# soil water pools
function neural_moisture(n::NODE, u0, ps, st, input, soildepth, perc, transp; dt = 1.0f0)
    
    st_model = Lux.StatefulLuxLayer{true}(n.model, ps, st)

    rhs(u, p, t) = perc - transp + st_model(vcat(u./soildepth, input), p) 
    
    prob = ODEProblem{false}(ODEFunction{false}(rhs), u0, n.tspan, ps)
    
    return solve(prob, n.solver; dt = dt)
end

function neural_moisture(n::NODE, u0, ps, st, input, soildepth, perc, transp, evapor; dt = 1.0f0)
    
    st_model = Lux.StatefulLuxLayer{true}(n.model, ps, st)

    rhs(u, p, t) = perc - transp - evapor + st_model(vcat(u./soildepth, input), p) 
    
    prob = ODEProblem{false}(ODEFunction{false}(rhs), u0, n.tspan, ps)
    
    return solve(prob, n.solver; dt = dt)
end

function neural_moisture(n::NODE, u0, ps, st, input, soildepth; dt = 1.0f0)
    
    st_model = Lux.StatefulLuxLayer{true}(n.model, ps, st)

    rhs(u, p, t) = st_model(vcat(u./soildepth, input), p) 
    
    prob = ODEProblem{false}(ODEFunction{false}(rhs), u0, n.tspan, ps)
    
    return solve(prob, n.solver; dt = dt)
end