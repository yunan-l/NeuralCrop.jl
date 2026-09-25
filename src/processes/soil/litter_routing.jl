"""
    litter_tillage!(soil, crop)

On sowing days, transfer the configured fraction of existing surface litter
to the incorporated litter pool. This follows LPJmL's `cultivate.c` ->
`tillage.c` order.
"""

function litter_tillage!(soil, crop)
    launch_custom!(
        litter_tillage_kernel!,
        soil_carbon_prognostic(soil).litter,
        size(soil_carbon_prognostic(soil).litter, 2),
        soil_nitrogen_prognostic(soil).litter,
        soil_management_input(soil).tillage_fraction,
        crop_events(crop).sowing,
        soil_management_fluxes(soil).tillage_carbon,
        soil_management_fluxes(soil).tillage_nitrogen,
    )
    return nothing
end

@kernel inbounds = true function litter_tillage_kernel!(
    carbon_litter::AbstractMatrix{T},
    nitrogen_litter::AbstractMatrix{T},
    tillage_fraction::AbstractMatrix{T},
    sowing_event::AbstractVector{S},
    tillage_carbon::AbstractVector{T},
    tillage_nitrogen::AbstractVector{T},
) where {T <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    event = sowing_event[cell]
    if event != 0
        carbon_1 = carbon_litter[1, cell]
        carbon_2 = carbon_litter[2, cell]
        carbon_3 = carbon_litter[3, cell]
        nitrogen_1 = nitrogen_litter[1, cell]
        nitrogen_2 = nitrogen_litter[2, cell]
        nitrogen_3 = nitrogen_litter[3, cell]
        for destination in 1:3
            carbon_routed = tillage_fraction[destination, 1] * carbon_1 +
                tillage_fraction[destination, 2] * carbon_2 +
                tillage_fraction[destination, 3] * carbon_3
            nitrogen_routed = tillage_fraction[destination, 1] * nitrogen_1 +
                tillage_fraction[destination, 2] * nitrogen_2 +
                tillage_fraction[destination, 3] * nitrogen_3
            carbon_litter[destination, cell] = carbon_routed
            nitrogen_litter[destination, cell] = nitrogen_routed
        end
        tillage_carbon[cell] = max(carbon_1 - carbon_litter[1, cell], zero(T))
        tillage_nitrogen[cell] = max(nitrogen_1 - nitrogen_litter[1, cell], zero(T))
    else
        tillage_carbon[cell] = zero(T)
        tillage_nitrogen[cell] = zero(T)
    end
end

"""
    litter_bioturbation!(soil; lpjmlparams=lpjmlparams)

Apply LPJmL's daily bioturbation transfer from surface (`agtop`) to
incorporated (`agsub`) litter. Carbon and nitrogen are moved together and the
operation is conservative for each cell.
"""

function litter_bioturbation!(soil;
                              lpjmlparams::LPJmLParams = lpjmlparams)
    launch_custom!(
        litter_bioturbation_kernel!,
        soil_carbon_prognostic(soil).litter,
        size(soil_carbon_prognostic(soil).litter, 2),
        soil_nitrogen_prognostic(soil).litter,
        soil_management_fluxes(soil).bioturbation_carbon,
        soil_management_fluxes(soil).bioturbation_nitrogen,
        kernel_constant(soil_carbon_prognostic(soil).litter, lpjmlparams),
    )
    return nothing
end

@kernel inbounds = true function litter_bioturbation_kernel!(
    carbon_litter::AbstractMatrix{T},
    nitrogen_litter::AbstractMatrix{T},
    bioturbation_carbon::AbstractVector{T},
    bioturbation_nitrogen::AbstractVector{T},
    fixed_parameters,
) where {T <: AbstractFloat}
    cell = @index(Global)
    lpjmlparams = kernel_value(fixed_parameters)
    fraction = T(lpjmlparams.bioturbate)
    carbon_flux = carbon_litter[SURFACE_LITTER, cell] * fraction
    nitrogen_flux = nitrogen_litter[SURFACE_LITTER, cell] * fraction
    bioturbation_carbon[cell] = carbon_flux
    bioturbation_nitrogen[cell] = nitrogen_flux
    carbon_litter[INCORPORATED_LITTER, cell] += carbon_flux
    carbon_litter[SURFACE_LITTER, cell] -= carbon_flux
    nitrogen_litter[INCORPORATED_LITTER, cell] += nitrogen_flux
    nitrogen_litter[SURFACE_LITTER, cell] -= nitrogen_flux
end

"""
Route today's harvested carbon residues through LPJmL's post-harvest tillage.

`harvest_crop.c` first adds shoot residues to `agtop` and roots to `bg`.
The harvested stand is then marked `KILL`; the same day's `killstand()` calls
`setaside()`, which calls `tillage()` when tillage is enabled. The root pool is
unchanged by the tillage matrix.
"""

function route_harvest_carbon_input!(soil, crop)
    launch_custom!(
        route_harvest_litter_kernel!,
        soil_carbon_prognostic(soil).litter,
        size(soil_carbon_prognostic(soil).litter, 2),
        soil_carbon_fluxes(soil).input,
        soil_management_input(soil).tillage_fraction,
        crop_events(crop).harvest,
        soil_management_fluxes(soil).tillage_carbon,
    )
    return nothing
end

"""Route today's harvested nitrogen residues through post-harvest tillage."""


function route_harvest_nitrogen_input!(soil, crop)
    launch_custom!(
        route_harvest_litter_kernel!,
        soil_nitrogen_prognostic(soil).litter,
        size(soil_nitrogen_prognostic(soil).litter, 2),
        soil_nitrogen_fluxes(soil).input,
        soil_management_input(soil).tillage_fraction,
        crop_events(crop).harvest,
        soil_management_fluxes(soil).tillage_nitrogen,
    )
    return nothing
end

@kernel inbounds = true function route_harvest_litter_kernel!(
    litter::AbstractMatrix{T},
    litter_input::AbstractMatrix{T},
    tillage_fraction::AbstractMatrix{T},
    harvest_event::AbstractVector{S},
    tillage_flux::AbstractVector{T},
) where {T <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    if harvest_event[cell] != 0
        # LPJmL routes the pools verbatim. In a negative-biomass failure the
        # mobile pool can be negative, and clipping it here would create carbon.
        litter_1 = litter[1, cell] + litter_input[1, cell]
        litter_2 = litter[2, cell] + litter_input[2, cell]
        litter_3 = litter[3, cell] + litter_input[3, cell]
        routed_1 = tillage_fraction[1, 1] * litter_1 +
            tillage_fraction[1, 2] * litter_2 + tillage_fraction[1, 3] * litter_3
        tillage_flux[cell] += max(litter_1 - routed_1, zero(T))
        litter[1, cell] = routed_1
        litter[2, cell] = tillage_fraction[2, 1] * litter_1 +
            tillage_fraction[2, 2] * litter_2 + tillage_fraction[2, 3] * litter_3
        litter[3, cell] = tillage_fraction[3, 1] * litter_1 +
            tillage_fraction[3, 2] * litter_2 + tillage_fraction[3, 3] * litter_3
    end
end
