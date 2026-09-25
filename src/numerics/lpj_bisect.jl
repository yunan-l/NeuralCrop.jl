"""
    lpj_bisect(f, xlow, xhigh; x_accuracy=0, y_accuracy=1e-3,
               max_iterations=30)

CPU reference implementation of LPJmL's bisection routine. If the stopping
criteria are not met, it returns the sampled point with the smallest absolute
function value, matching LPJmL rather than requiring a strictly bracketed root.

Returns `(root, iterations)`.
"""
function lpj_bisect(f,
                    xlow::T,
                    xhigh::T;
                    x_accuracy::T = zero(T),
                    y_accuracy::T = T(1e-3),
                    max_iterations::Integer = 30) where {T <: AbstractFloat}
    ylow = f(xlow)
    xmin = (xlow + xhigh) * T(0.5)
    ymin = typemax(T)

    for iteration in 0:(max_iterations - 1)
        xmid = (xlow + xhigh) * T(0.5)
        if xhigh - xlow < x_accuracy
            return xmid, iteration
        end

        ymid = f(xmid)
        if abs(ymid) < ymin
            ymin = abs(ymid)
            xmin = xmid
        end
        if abs(ymid) < y_accuracy
            return xmid, iteration
        end

        if ylow * ymid <= zero(T)
            xhigh = xmid
        else
            xlow = xmid
            ylow = ymid
        end
    end

    return xmin, max_iterations
end
