function petpar!(pet::PetPar,
                 day::Int64,
                 lat::AbstractArray{T},
                 temp::AbstractArray{T},
                 lwnet::AbstractArray{T},
                 swdown::AbstractArray{T};
                 dayseconds = 86400
) where {T <: AbstractFloat}

    
    delta = Float32(deg2rad(-23.4 * cos(2 * π * (day + 10) / 365)))
    u = Float32.(sin.(deg2rad.(lat)) * sin(delta))
    v = Float32.(cos.(deg2rad.(lat)) * cos(delta))
    
    backend = KernelAbstractions.get_backend(lat)
    
    kernel = daylength_kernel!(backend)
    
    kernel(pet.daylength, u, v, ndrange=length(lat))
    
    KernelAbstractions.synchronize(backend)
    
    swnet = (1 .- pet.albedo) .* swdown
    
    pet.par .= dayseconds .* swdown ./ 2
    
    s = 2.503f6 * exp.(17.269f0 * temp ./ (237.3f0 .+ temp)) ./ ((237.3f0 .+ temp).^2)

    gamma_t = 65.05f0 .+ 0.064f0 * temp
    lambda = 2.495f6 .- 2380f0 * temp

    pet.eeq .= dayseconds * (s ./ (s .+ gamma_t) ./ lambda) .* (swnet .+ lwnet .* (pet.daylength / 24))
    
    # idx = pet.eeq .< 0
    # pet.eeq[idx] .= zero(T)
    # pet.eeq .= ifelse.(pet.eeq .< 0, zero(T), pet.eeq)
    pet.eeq .= max.(pet.eeq, zero(T))

    ## check equilibrium evapotranspiration
    pet.eeq .= min.(pet.eeq, 15.0f0) ##  set an upper bound for pet.eeq to avoid extreme values to stop GPU computing

end

@kernel function daylength_kernel!(pet_daylength::AbstractArray{T},
                                   u::AbstractArray{T},
                                   v::AbstractArray{T}
) where {T <: AbstractFloat}
    
    cell = @index(Global)
    
    if u[cell] >= v[cell]
        pet_daylength[cell] = 24
    elseif u[cell] <= -v[cell]
        pet_daylength[cell] = 0
    else
        hh = acos(-u[cell] / v[cell])
        pet_daylength[cell] = 24 * hh * (1 / π)
    end
end

# for one cft
function apar_crop!(PFT::PftParameters,
                    crop::Crop,
                    pet::PetPar
)
    
    @unpack name, lightextcoeff, albedo_leaf, alphaa  = PFT
    
    # crop.fpar .= (1 .- exp.(-lightextcoeff * max.(0.0f0, crop.lai .- crop.lai_nppdeficit))) # if maize, crop.fpar = min.(1.0f0, max.(0.0f0, 0.2558f0 * max.(0.01f0, crop.lai .- crop.lai_nppdeficit) .- 0.0024f0))
    crop.fpar .= 1 .- exp.(-lightextcoeff * max.(0.0f0, crop.lai))
    
    crop.apar .= pet.par * (1 - albedo_leaf) * alphaa .* crop.fpar

end 

function apar_crop_maize!(PFT::PftParameters,
                          crop::Crop,
                          pet::PetPar
)
    
    @unpack name, lightextcoeff, albedo_leaf, alphaa  = PFT
    
    # crop.fpar = min.(1.0f0, max.(0.0f0, 0.2558f0 * max.(0.01f0, crop.lai .- crop.lai_nppdeficit) .- 0.0024f0))
    crop.fpar = min.(1.0f0, max.(0.0f0, 0.2558f0 * max.(0.01f0, crop.lai) .- 0.0024f0))
    
    crop.apar .= pet.par * (1 - albedo_leaf) * alphaa .* crop.fpar

end 