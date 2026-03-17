function deg2rad(x) # Convert degree -> radian
    return x * π / 180.0
end
      
    
function ppm2Pa(co2::AbstractArray{T}) where {T <: AbstractFloat} # Convert ppmv --> Pa
        
    return co2 * T(1e-1)
end

function ppm2bar(co2::AbstractArray{T}) where {T <: AbstractFloat} # Convert ppmv --> bar
        
    return co2 * T(1e-6)
end
    
function hour2day(hour::AbstractArray{T}) where {T <: AbstractFloat} # Convert hour --> day
        
    return hour * T(0.04166666666666666666)
end

function hour2sec(hour::AbstractArray{T}) where {T <: AbstractFloat} # Convert hour --> sec
        
    return hour * T(3600)
end
    
function degCtoK(deg::AbstractArray{T}) where {T <: AbstractFloat} # deg C --> Kelvin
        
    return deg .+ T(273.15)
end