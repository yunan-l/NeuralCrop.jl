"""Convert degrees to radians."""
function deg2rad(x) # Convert degree -> radian
    return x * π / 180.0
end
      
    
"""Convert CO₂ concentration from ppmv to Pa."""
function ppm2Pa(co2::AbstractArray{T}) where {T <: AbstractFloat} # Convert ppmv --> Pa
        
    return co2 * T(1e-1)
end

"""Convert CO₂ concentration from ppmv to bar."""
function ppm2bar(co2::AbstractArray{T}) where {T <: AbstractFloat} # Convert ppmv --> bar
        
    return co2 * T(1e-6)
end
    
"""Convert hours to day fraction."""
@inline function hour2day(hour::T) where {T <: AbstractFloat}
    return hour * T(0.04166666666666666666)
end

function hour2day(hour::AbstractArray{T}) where {T <: AbstractFloat} # Convert hour --> day
    return hour * T(0.04166666666666666666)
end

"""Convert hours to seconds."""
@inline function hour2sec(hour::T) where {T <: AbstractFloat}
    return hour * T(3600)
end

function hour2sec(hour::AbstractArray{T}) where {T <: AbstractFloat} # Convert hour --> sec
        
    return hour * T(3600)
end
    
"""Convert temperature from degrees Celsius to Kelvin."""
function degCtoK(deg::AbstractArray{T}) where {T <: AbstractFloat} # deg C --> Kelvin
        
    return deg .+ T(273.15)
end
