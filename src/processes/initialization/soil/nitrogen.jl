"""
Mineral and organic soil nitrogen pools, inputs, and decomposition fluxes.

The shared fixed normalized post-spin-up routing distributions are stored in
`soil.decomposition`; the `litter_to_*` arrays store the actual daily
layer-resolved retained-N fluxes.
"""
mutable struct SoilNitrogen{A, L, M}
    nitrate::M                      # Soil nitrate stock by layer (gN m⁻²).
    ammonium::M                     # Soil ammonium stock by layer (gN m⁻²).
    litter::L                       # Organic nitrogen in the three litter classes (gN m⁻²).
    decomposed_litter::L            # Litter nitrogen decomposed today (gN m⁻² day⁻¹).
    fast::M                         # Fast soil-organic-nitrogen stock by layer (gN m⁻²).
    slow::M                         # Slow soil-organic-nitrogen stock by layer (gN m⁻²).
    decomposed_fast::M              # Fast organic-N decomposed today (gN m⁻² day⁻¹).
    decomposed_slow::M              # Slow organic-N decomposed today (gN m⁻² day⁻¹).
    input::L                        # Current-day N input to litter classes (gN m⁻² day⁻¹).
    litter_to_fast::M               # Retained litter N routed to fast pool today (gN m⁻² day⁻¹).
    litter_to_slow::M               # Retained litter N routed to slow pool today (gN m⁻² day⁻¹).
    litter_response::A              # Environmental decomposition multiplier by litter class (0–1+).
    mineralization::M               # Net organic-N mineralized to NH₄ by layer (gN m⁻² day⁻¹).
    immobilization::M               # Mineral N immobilized into organic pools by layer (gN m⁻² day⁻¹).
    nitrification::M                # NH₄ converted by nitrification by layer (gN m⁻² day⁻¹).
    n2o_nitrification::M            # N₂O-N emitted during nitrification by layer (gN m⁻² day⁻¹).
    denitrification::M              # NO₃ consumed by denitrification by layer (gN m⁻² day⁻¹).
    n2o_denitrification::M          # N₂O-N emitted during denitrification by layer (gN m⁻² day⁻¹).
    n2_denitrification::M           # N₂-N emitted during denitrification by layer (gN m⁻² day⁻¹).
    volatilization::A               # NH₃-N volatilized from the soil surface (gN m⁻² day⁻¹).
    leaching::A                     # Mineral N removed by bottom drainage (gN m⁻² day⁻¹).
end

init_soil_nitrogen(cell_size::Int, device; kwargs...) =
    init_soil_nitrogen(Float32, cell_size, device; kwargs...)
function init_soil_nitrogen(::Type{T}, cell_size::Int, device;
                            litter_layers::Int = 3,
                            soil_layers::Int = 5) where {T <: AbstractFloat}
    litter_state() = device(zeros(T, litter_layers, cell_size))
    layer_state() = device(zeros(T, soil_layers, cell_size))
    return SoilNitrogen(
        layer_state(), layer_state(),
        litter_state(), litter_state(),
        layer_state(), layer_state(), layer_state(), layer_state(),
        litter_state(), layer_state(), layer_state(),
        device(zeros(T, litter_layers)),
        layer_state(), layer_state(), layer_state(), layer_state(),
        layer_state(), layer_state(), layer_state(),
        device(zeros(T, cell_size)),
        device(zeros(T, cell_size)),
    )
end
