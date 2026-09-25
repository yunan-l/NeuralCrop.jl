"""Persistent seasonal crop-water memory."""
mutable struct CropWaterState{A}
    demand_sum::A # Seasonal accumulated transpiration demand (mm).
    supply_sum::A # Seasonal accumulated transpiration supply (mm).
    sufficiency::A # Prior-day water-sufficiency multiplier needed by next-day canopy growth (0–1).
end

"""Current-day crop-water fluxes."""
mutable struct CropWaterFluxes{A, M}
    interception::A       # Rainfall intercepted and evaporated by the canopy (mm day⁻¹).
    transpiration_layer::M # Root-water uptake/transpiration from each soil layer (mm day⁻¹).
end

function init_crop_water_state(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    float_state() = device(zeros(T, cell_size))
    return CropWaterState(float_state(), float_state(), device(ones(T, cell_size)))
end

function init_crop_water_fluxes(::Type{T}, cell_size::Int, device;
                                soil_layers::Int = 5) where {T <: AbstractFloat}
    float_flux() = device(zeros(T, cell_size))
    return CropWaterFluxes(
        float_flux(),
        device(zeros(T, soil_layers, cell_size)),
    )
end
