"""
z_score_norm(x)

Apply z-score normalization along feature dimensions.
Returns normalized data, mean, and standard deviation.
"""
function z_score_norm(x::AbstractArray{T}; by_columns::Bool = false, precipitation::Bool = false) where {T <: AbstractFloat}
    
    if by_columns
        μ = vec(mean(x; dims = 1))
        σ = vec(std(x; dims = 1))
        x_norm = similar(x)

        for i in axes(x, 2)
            if σ[i] == 0
                x_norm[:, i] .= 0.0f0
            else
                x_norm[:, i] = (x[:, i] .- μ[i]) ./ σ[i]
            end
        end
    else
        if precipitation
            x = log1p.(x)
        end
        μ = mean(x)
        σ = std(x)
        if σ == 0
            x_norm = zeros(size(x))
        else
            x_norm = (x .- μ) ./ σ
        end
    end

    return x_norm, μ, σ
end

"""
apply_z_score(x, μ, σ)

Normalize `x` using precomputed mean `μ` and standard deviation `σ`.
"""
function apply_z_score(x::AbstractArray{T}, μ, σ; by_columns::Bool = false, precipitation::Bool = false) where {T <: AbstractFloat}
    
    if by_columns
        x_norm = similar(x)
        for i in axes(x, 2)
            if σ[i] == 0
                x_norm[:, i] .= 0.0f0
            else
                x_norm[:, i] = (x[:, i] .- μ[i]) ./ σ[i]
            end
        end
    else
        if precipitation
            x = log1p.(x)
        end
        if σ == 0
            x_norm = zeros(size(x))
        else
            x_norm = (x .- μ) ./ σ
        end
    end
    return x_norm
end


function min_max_one_dimension(x)
    x_min = minimum(x; dims = 1)
    x_max = maximum(x; dims = 1)
    if x_max[1] == x_min[1]
        x_norm .= 0.0f0
        x_max[1] = x_min[1] + 1.0f0
    else
        x_norm = (x .- x_min) ./ (x_max .- x_min)
    end

    return x_norm, x_min, x_max
end

"""
min_max_norm(x)

Apply min-max normalization and return normalized data plus min/max statistics.
"""
function min_max_norm(x::AbstractArray{T}; by_columns::Bool = false) where {T <: AbstractFloat}
    
    if by_columns
        x_min = vec(minimum(x; dims=1))
        x_max = vec(maximum(x; dims=1))
        x_norm = similar(x)

        for i in axes(x, 2)
            if x_max[i] == x_min[i]
                x_norm[:, i] .= 0.0f0
            elseif abs(x_max[i]) < 1f-3
                x_norm[:, i] = x[:, i]
                x_max[i] = 1.0f0
                x_min[i] = 0.0f0
            elseif (x_max[i] - x_min[i]) < 1.0
                x_min[i] = 0.0f0
                x_norm[:, i] = (x[:, i] .- x_min[i]) ./ (x_max[i] - x_min[i])
            else
                x_norm[:, i] = (x[:, i] .- x_min[i]) ./ (x_max[i] - x_min[i])
            end
        end
    else
        x_min = minimum(x)
        x_max = maximum(x)
        if x_max == x_min
            x_norm = zeros(size(x))
        elseif abs(x_max) < 1f-3
            x_norm = x
            x_max = 1.0f0
            x_min = 0.0f0
        elseif (x_max - x_min) < 1.0
            x_min = 0.0f0
            x_norm = (x .- x_min) ./ (x_max - x_min)
        else
            x_norm = (x .- x_min) ./ (x_max - x_min)
        end
    end
    
    return x_norm, x_min, x_max
end
