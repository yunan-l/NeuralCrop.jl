"""
    nitrogen_deposition!(soil, nitrate_deposition, ammonium_deposition)

Add externally supplied daily atmospheric nitrate- and ammonium-N deposition
to the upper soil layer. Inputs use `gN m⁻² day⁻¹`, matching LPJmL's daily
climate fields. The caller supplies zeros when deposition forcing is absent.
"""
function nitrogen_deposition!(state::ModelState,
                               nitrate_deposition::AbstractVector{T},
                               ammonium_deposition::AbstractVector{T}) where {T <: AbstractFloat}
    mineral_nitrogen = soil_nitrogen_prognostic(state)
    return nitrogen_deposition!(
        mineral_nitrogen.nitrate, mineral_nitrogen.ammonium,
        nitrate_deposition, ammonium_deposition,
    )
end

function nitrogen_deposition!(soil::Soil,
                               nitrate_deposition::AbstractVector{T},
                               ammonium_deposition::AbstractVector{T}) where {T <: AbstractFloat}
    return nitrogen_deposition!(
        soil.nitrogen.nitrate, soil.nitrogen.ammonium,
        nitrate_deposition, ammonium_deposition,
    )
end

function nitrogen_deposition!(nitrate::AbstractMatrix{T},
                               ammonium::AbstractMatrix{T},
                               nitrate_deposition::AbstractVector{T},
                               ammonium_deposition::AbstractVector{T}) where {T <: AbstractFloat}
    launch_custom!(
        nitrogen_deposition_kernel!,
        nitrate,
        size(nitrate, 2),
        ammonium,
        nitrate_deposition,
        ammonium_deposition,
    )
    return nothing
end

@kernel inbounds = true function nitrogen_deposition_kernel!(
    nitrate::AbstractMatrix{T},
    ammonium::AbstractMatrix{T},
    nitrate_deposition::AbstractVector{T},
    ammonium_deposition::AbstractVector{T},
) where {T <: AbstractFloat}
    cell = @index(Global)
    ammonium[1, cell] += ammonium_deposition[cell]
    nitrate[1, cell] += nitrate_deposition[cell]
end
