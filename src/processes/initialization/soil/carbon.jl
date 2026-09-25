"""
Soil carbon inputs, pools, decomposition fluxes, and respiration.

The shared fixed post-spin-up vertical routing distributions are stored in
`soil.decomposition`. `litter_to_fast` and `litter_to_slow` are the daily
layer-resolved fluxes after applying fastfrac and the retained carbon fraction.
"""
mutable struct SoilCarbon{A, L, M}
    input::L                       # Current-day C input to surface/incorporated/root litter classes (gC m⁻² day⁻¹).
    litter::L                      # Carbon stock in the three litter classes (gC m⁻²).
    decomposed_litter::L           # Litter carbon decomposed today by class (gC m⁻² day⁻¹).
    fast::M                        # Fast soil-organic-carbon stock by soil layer (gC m⁻²).
    slow::M                        # Slow soil-organic-carbon stock by soil layer (gC m⁻²).
    decomposed_fast::M             # Fast-pool carbon decomposed today (gC m⁻² day⁻¹).
    decomposed_slow::M             # Slow-pool carbon decomposed today (gC m⁻² day⁻¹).
    litter_to_fast::M              # Retained litter carbon routed to fast pool today (gC m⁻² day⁻¹).
    litter_to_slow::M              # Retained litter carbon routed to slow pool today (gC m⁻² day⁻¹).
    litter_response::A             # Environmental decomposition multiplier for each litter class (0–1+).
    heterotrophic_respiration::A   # Total litter plus soil heterotrophic respiration (gC m⁻² day⁻¹).
end

init_soil_carbon(cell_size::Int, device; kwargs...) =
    init_soil_carbon(Float32, cell_size, device; kwargs...)
function init_soil_carbon(::Type{T}, cell_size::Int, device;
                          litter_layers::Int = 3,
                          soil_layers::Int = 5) where {T <: AbstractFloat}
    litter_state() = device(zeros(T, litter_layers, cell_size))
    layer_state() = device(zeros(T, soil_layers, cell_size))
    return SoilCarbon(
        litter_state(), litter_state(), litter_state(),
        layer_state(), layer_state(), layer_state(), layer_state(),
        layer_state(), layer_state(),
        device(zeros(T, litter_layers)),
        device(zeros(T, cell_size)),
    )
end
