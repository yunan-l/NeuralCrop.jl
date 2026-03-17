function temp_stress(PFT::PftParameters,
                     photopar::PhotoPar,
                     pet::PetPar,
                     photos::Photos,
                     k::K,
                     temp::AbstractArray{T}
) where {T <: AbstractFloat}

    backend = get_backend(temp)
    
    kernel = temp_stress_kernel!(backend)
    
    kernel(pet.daylength, temp, photos.tstress, k, photopar, PFT, ndrange=length(photos.tstress))
    
    synchronize(backend)
  
end


@kernel function temp_stress_kernel!(pet_daylength::AbstractArray{T},           
                                     temp::AbstractArray{T},           
                                     photos_tstress::AbstractArray{T},            
                                     k::K,
                                     photopar::PhotoPar,
                                     PFT::PftParameters
) where {T <: AbstractFloat}
    
    cell = @index(Global)

    @unpack path, temp_co2, temp_photos = PFT
    @unpack k1, k2, k3 = k
    @unpack tmc3, tmc4 = photopar
    
    if pet_daylength[cell] < 0.01 || (path == 1 && temp[cell] > tmc3) || (path == 2 && temp[cell] > tmc4) # path == 1 : C3; path == 2 : C4
        photos_tstress[cell] = zero(T)
    else
        if temp[cell] < temp_co2.high
            low = 1 / (1 + exp(k1 * (k2 - temp[cell])))
            high = 1 - 0.01 .* exp(k3 * (temp[cell] - temp_photos.high))
            photos_tstress[cell] = T(low * high)
        else
            photos_tstress[cell] = zero(T)
        end
    end
end