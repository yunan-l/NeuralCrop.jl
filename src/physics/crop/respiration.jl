function respiration!(crop::Crop,
                      PFT::PftParameters,
                      param::LPJmLParam,
                      temp::AbstractArray{T},
                      assim::AbstractArray{T}
) where {T <: AbstractFloat}
    
    @unpack respcoeff, nc_ratio = PFT
    @unpack k, r_growth, e0, temp_response = param
    
    # kernel based
    gtemp_air = similar(temp)
    backend = get_backend(temp)
    kernel = temp_response_kernel!(backend)
    kernel(temp_response, e0, temp, gtemp_air, ndrange=length(assim))
    synchronize(backend)
    
    # unlimited nitrogen
    rosoresp = crop.rootc * respcoeff * k * nc_ratio.root .* gtemp_air .+ crop.stoc * respcoeff * k * nc_ratio.sto .* gtemp_air
    presp = crop.poolc * respcoeff * k * nc_ratio.pool .* gtemp_air
    gresp = (assim .- rosoresp .- presp) * r_growth

    # # differentiation based
    # # gate = max.(temp .+ T(40.0), T(0.0)) ./ (max.(temp .+ T(40.0), T(1e-5)))
    # gate = sigmoid.(T(10.0) * (temp .+ T(40.0)))
    # gtemp_air = gate .* exp.(e0 * (one(T) / (temp_response + T(10.0)) .- one(T) ./ (temp .+ temp_response)))
    # rosoresp = crop.rootc * respcoeff * k .* (crop.rootn ./ (crop.rootc .+ T(1e-5))) .* gtemp_air .+ crop.stoc * respcoeff * k .* (crop.ston ./ (crop.stoc .+ T(1e-5))) .* gtemp_air
    # presp = crop.poolc * respcoeff * k .* (crop.pooln ./ (crop.poolc .+ T(1e-5))) .* gtemp_air
    # gresp = (assim .- rosoresp .- presp) * r_growth
    # gresp = ifelse.(gresp .< zero(T), zero(T), gresp)

    
    crop.resp = (rosoresp .+ presp .+ gresp) .* crop.isgrowing

end
    
    
@kernel function temp_response_kernel!(temp_response::T,
                                       e0::T,
                                       temp::AbstractArray{T},
                                       gtemp_response::AbstractArray{T}
                                       
            
) where {T <: AbstractFloat}
    
    cell = @index(Global)
    
    if temp[cell] >= -40.0
        gtemp_response[cell] = exp(e0 * (one(T) / (temp_response + T(10.0)) - one(T) / (temp[cell] + temp_response)))
    else
        gtemp_response[cell] = zero(T)
    end
end