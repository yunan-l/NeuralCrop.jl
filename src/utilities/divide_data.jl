# only used for one grid cell running
function divide_data_one_dimension(data::AbstractArray{T}, 
                                   sdate::AbstractArray{S}, 
                                   hdate::AbstractArray{S}
) where {T <: AbstractFloat, S <: Integer}

    out = []
    for i in 1:length(sdate)
        data_ = data[sdate[i]:hdate[i]-1]
        push!(out, data_)
    end
    return out
end

function divide_data_dimensions(data::AbstractArray{T}, 
                                sdate::AbstractArray{S}, 
                                hdate::AbstractArray{S}
) where {T <: AbstractFloat, S <: Integer}

    out = []
    for i in 1:length(sdate)
        data_ = data[:, sdate[i]:hdate[i]-1]
        push!(out, data_)
    end
    return out
end

function divide_data_one_dimension(data::AbstractArray{T}, 
                                   sdate::AbstractArray{S}
) where {T <: AbstractFloat, S <: Integer}

    out = []
    for i in 2:length(sdate)
        data_ = data[sdate[i-1]-1:sdate[i]-2]
        push!(out, data_)
    end
    return out
end

function divide_data_dimensions(data::AbstractArray{T}, 
                                sdate::AbstractArray{S}
) where {T <: AbstractFloat, S <: Integer}

    out = []
    for i in 2:length(sdate)
        data_ = data[:, sdate[i-1]-1:sdate[i]-2]
        push!(out, data_)
    end
    return out
end