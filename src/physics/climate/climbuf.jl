using CUDA

function annual_climbuf!(daily_temp::AbstractArray{T},
                         climbuf::ClimBuf,
                         PFT::PftParameters,
                         device;
                         n = 5,
                         kk = 0.05
) where {T <: AbstractFloat}
    # Calculate the average temperature for each month.
    # update_monthly
    # length(daily_temp) = 365
    # n = 5, the first n coldest months
    # kk is to rescale the 20-year average monthly temprerature
    
    monthlytemp!(daily_temp, climbuf.mtemp, device)
    
    backend = KernelAbstractions.get_backend(daily_temp)

    kernel = climbuf_mtemp20_kernel!(backend)

    kernel(climbuf.mtemp20, climbuf.mtemp, kk, ndrange=(size(climbuf.mtemp20, 1), size(climbuf.mtemp20, 2)))

    KernelAbstractions.synchronize(backend)
    # climbuf.mtemp20 .= ifelse.(climbuf.mtemp20 .< -9998, climbuf.mtemp, (1 - kk) * climbuf.mtemp20 .+ kk * climbuf.mtemp)
    
    # getmintemp20_n!(climbuf, n)
    climbuf.min_temp .= sort(climbuf.mtemp20, dims=1)[1:n, :] # Array to store n coldest months
    
    backend = KernelAbstractions.get_backend(daily_temp)

    kernel = climbuf_V_req_a_kernel!(backend)

    kernel(climbuf.min_temp, climbuf.V_req_a, PFT, n, ndrange=length(climbuf.V_req_a))

    KernelAbstractions.synchronize(backend)
    
    # for m = 1:n
    #     if climbuf.min_temp[m] <= PFT.tv_opt.low && climbuf.min_temp[m]> -9999
    #         climbuf.V_req_a += PFT.pvd_max/ n
    #     elseif climbuf.min_temp[m] > PFT.tv_opt.low && climbuf.min_temp[m] < PFT.tv_opt.high
    #         climbuf.V_req_a += PFT.pvd_max / n * (1-(climbuf.min_temp[m] - PFT.tv_opt.low) / (PFT.tv_opt.high - PFT.tv_opt.low))
    #     end
    # end
    
    backend = KernelAbstractions.get_backend(daily_temp)

    kernel = climbuf_V_req_kernel!(backend)

    kernel(climbuf.V_req, climbuf.V_req_a, kk, ndrange=length(climbuf.V_req))
    
    KernelAbstractions.synchronize(backend)
    # climbuf.V_req .= ifelse.(climbuf.V_req .< -9998, climbuf.V_req_a, (1 - kk) * climbuf.V_req .+ kk .* climbuf.V_req_a)

    climbuf.atemp_mean = vec(mean(daily_temp, dims = 1))

end


@kernel function climbuf_mtemp20_kernel!(climbuf_mtemp20::AbstractArray{T},
                                         climbuf_mtemp::AbstractArray{T},
                                         kk
) where {T <: AbstractFloat}
    
    month, cell = @index(Global, NTuple)
    
    if climbuf_mtemp20[month, cell] < -9998
        climbuf_mtemp20[month, cell] = climbuf_mtemp[month, cell]
    else
        climbuf_mtemp20[month, cell] = (1 - kk) * climbuf_mtemp20[month, cell] + kk * climbuf_mtemp[month, cell]
    end
    
end


@kernel function climbuf_V_req_a_kernel!(climbuf_min_temp::AbstractArray{T},
                                         climbuf_V_req_a::AbstractArray{T},
                                         PFT::PftParameters,
                                         n
) where {T <: AbstractFloat}
    
    cell = @index(Global)

    @unpack tv_opt, pvd_max = PFT
    
    sum_v_req = zero(T)

    for i in 1:n
        if climbuf_min_temp[i, cell] <= tv_opt.low && climbuf_min_temp[i, cell]> -9999
            sum_v_req += pvd_max / n
        elseif climbuf_min_temp[i, cell] > tv_opt.low && climbuf_min_temp[i, cell] < tv_opt.high
            sum_v_req += pvd_max / n * (1 - (climbuf_min_temp[i, cell] - tv_opt.low) / (tv_opt.high - tv_opt.low))
        end
    end
    
    climbuf_V_req_a[cell] = sum_v_req

end


@kernel function climbuf_V_req_kernel!(climbuf_V_req::AbstractArray{T},
                                       climbuf_V_req_a::AbstractArray{T},
                                       kk
) where {T <: AbstractFloat}
    
    cell = @index(Global)
    
    if climbuf_V_req[cell] < -9998
        climbuf_V_req[cell] = climbuf_V_req_a[cell]
    else
        climbuf_V_req[cell] = (1 - kk) * climbuf_V_req[cell] + kk * climbuf_V_req_a[cell]
    end
    
end

# function getmintemp20_n!(climbuf::ClimBuf,
#                          n::Int
# )
#     """
#     Calculates the n coldest months from the climate buffer and returns their values.

#     Args:
#         climbuf: A dictionary containing climate data, specifically `:mtemp20` for monthly temperatures.
#         n: The number of coldest months to extract.

#     Return:
#         A vector containing the n coldest monthly temperatures.
#     """
#     climbuf.min_temp = sort(climbuf.mtemp20, dims=1)[1:n, :] # Array to store n coldest months
    
#     # for i in 1:n
#     #     index = argmin(temp[i:NMONTH]) + (i - 1)
#     #     min_temp[i] = temp[index]
#     #     # Swap the values
#     #     temp[i], temp[index] = temp[index], temp[i]
#     # end
# end

function monthlytemp!(daily_temp::AbstractArray{T},
                      climbuf_mtemp::AbstractArray{T},
                      device
) where {T <: AbstractFloat}
    """
    Calculate the average temperature for each month.

    Args:
        daily_temps::Vector{Float64}: Daily temperature data (length is 365)

    Return:
        A vector of length 12 representing the average temperature for each month.
    """
    ndaymonth = device([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]) #|> gdev
    start_indices = device([1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]) #|> gdev
    
    # start_indices = cumsum(vcat(1, ndaymonth[1:end-1]))
    # mtemp = similar(daily_temp, (12, cell_size))  # Store mean temperatures for each month
#     start_idx = 1  # Start index for each month in daily_temps

#     for month = 1:12
#         end_idx = start_idx + ndaymonth[month] - 1  # End index for the current month
#         mtemp[month] = mean(daily_temp[start_idx : end_idx])  # Calculate the monthly average
#         start_idx = end_idx + 1  # Update start index for the next month
#     end
    
    backend = KernelAbstractions.get_backend(daily_temp)
    
    kernel = monthlytemp_kernel!(backend)
    
    kernel(daily_temp, climbuf_mtemp, ndaymonth, start_indices, ndrange=(size(ndaymonth, 1), size(climbuf_mtemp, 2)))
    
    KernelAbstractions.synchronize(backend)
    
end


@kernel function monthlytemp_kernel!(daily_temp::AbstractArray{T},
                                     climbuf_mtemp::AbstractArray{T}, 
                                     ndaymonth::AbstractArray{S},  
                                     start_indices::AbstractArray{S}
) where {T <: AbstractFloat, S <: Integer}
    
    month, cell = @index(Global, NTuple)
    start_idx = start_indices[month]
    days = ndaymonth[month]

    sum_temp = zero(T)
    
    for i in 0:(days - 1)
        sum_temp += daily_temp[start_idx + i, cell]
    end

    climbuf_mtemp[month, cell] = sum_temp / days 

end


function daily_climbuf!(temp::AbstractArray{T},
                        climbuf_temp::AbstractArray{T}
) where {T <: AbstractFloat}

    backend = KernelAbstractions.get_backend(temp)

    kernel = daily_climbuf_kernel!(backend)

    kernel(temp, climbuf_temp, ndrange=length(temp))

    KernelAbstractions.synchronize(backend)

end


@kernel function daily_climbuf_kernel!(temp::AbstractArray{T},
                                       climbuf_temp::AbstractArray{T};
                                       NDAYS = 31
) where {T <: AbstractFloat}

    cell = @index(Global)

    for day in 2:NDAYS
        climbuf_temp[day-1, cell] = climbuf_temp[day, cell]
    end
    climbuf_temp[NDAYS, cell] = temp[cell]

end