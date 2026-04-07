struct SciMLEuler
end 

function solve(prob::SciMLBase.AbstractDEProblem, solver::SciMLEuler; kwargs...)

    u = prob.u0
    f = prob.f
    p = prob.p
    t = prob.tspan[1]
    
    dt = kwargs[:dt]
    
    u = @muladd u + dt .* f(u, p, t)

    return u
end

struct SciMLEuler_litc
end 

function solve(prob::SciMLBase.AbstractDEProblem, solver::SciMLEuler_litc; kwargs...)

    u = prob.u0
    f = prob.f
    p = prob.p
    t = prob.tspan[1]
    
    dt = kwargs[:dt]
    
    u = @muladd u + dt .* f(u, p, t)

    return u, -f(u, p, t)
end


struct SciMLEuler_soilc
end 

function solve(prob::SciMLBase.AbstractDEProblem, solver::SciMLEuler_soilc; kwargs...)

    u = prob.u0
    f = prob.f
    p = prob.p
    t = prob.tspan[1]
    
    A_trans = kwargs[:A_trans]
    c_input = kwargs[:c_input]
    dt = kwargs[:dt]
    
    u =  @muladd u + dt .* (A_trans .* c_input + f(u, p, t))

    return u, -f(u, p, t)
end