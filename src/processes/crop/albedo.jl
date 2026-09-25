"""
    albedo!(CFT, crop, soil, pet; maize=false)

Compute LPJmL-style effective surface albedo from green canopy absorption,
surface-litter cover, bare soil, and the snow state present at the start of the
day. Surface-litter cover is rebuilt directly from its carbon stock so event
routing and restart initialization cannot leave radiation on a stale cache.
"""

function albedo!(CFT::CFTParameters,
                 crop,
                 soil,
                 pet::PetPar;
                 maize::Bool = false,
                 soil_albedo = 0.3f0,
                 snow_albedo = 0.65f0,
                 litter_carbon_fraction = 0.42f0)
    launch_1D!(
        albedo_kernel!,
        pet.albedo,
        crop_canopy_auxiliary(crop).albedo,
        crop_prognostic(crop).canopy.lai,
        crop_prognostic(crop).canopy.lai_npp_deficit,
        crop_prognostic(crop).phenology.is_growing,
        soil_carbon_prognostic(soil).litter,
        soil_snow_prognostic(soil).height,
        soil_snow_prognostic(soil).fraction,
        kernel_constant(
            pet.albedo,
            (; CFT, soil_albedo, snow_albedo, litter_carbon_fraction),
        ),
        maize,
    )
    return nothing
end

@kernel inbounds = true function albedo_kernel!(
    pet_albedo::AbstractVector{T},
    canopy_albedo::AbstractVector{T},
    lai::AbstractVector{T},
    lai_npp_deficit::AbstractVector{T},
    is_growing::AbstractVector{S},
    carbon_litter::AbstractMatrix{T},
    snow_height::AbstractVector{T},
    snow_fraction::AbstractVector{T},
    fixed_parameters,
    maize::Bool,
) where {T <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    surface_parameters = kernel_value(fixed_parameters)
    CFT = surface_parameters.CFT
    light_extinction = T(CFT.lightextcoeff)
    leaf_albedo = T(CFT.albedo_leaf)
    litter_albedo = T(CFT.albedo_litter)
    fpc = T(CFT.fpc)
    soil_albedo = T(surface_parameters.soil_albedo)
    snow_albedo = T(surface_parameters.snow_albedo)
    litter_carbon_fraction = T(surface_parameters.litter_carbon_fraction)
    actual_lai = max(zero(T), lai[cell] - lai_npp_deficit[cell])
    green_fraction = maize ?
        clamp(T(0.2558) * max(T(0.01), actual_lai) - T(0.0024), zero(T), one(T)) :
        one(T) - exp(-light_extinction * actual_lai)
    litter_dry_matter = max(carbon_litter[1, cell], zero(T)) /
        litter_carbon_fraction
    litter_cover = one(T) - exp(-T(6e-3) * litter_dry_matter)
    background_fraction = one(T) - green_fraction
    snow_present = snow_height[cell] > zero(T)

    green_component = green_fraction *
        (snow_present ? snow_albedo : leaf_albedo)
    litter_component = litter_cover * background_fraction *
        (snow_present ? snow_albedo : litter_albedo)
    soil_background = (one(T) - litter_cover) * background_fraction
    soil_component = snow_present ?
        soil_background * snow_fraction[cell] * snow_albedo :
        soil_background * soil_albedo
    crop_surface = green_component + litter_component + soil_component
    bare_surface = snow_fraction[cell] * snow_albedo +
        (one(T) - snow_fraction[cell]) * soil_albedo

    if is_growing[cell] != zero(S)
        canopy_albedo[cell] = crop_surface
        pet_albedo[cell] = crop_surface + max(one(T) - fpc, zero(T)) * bare_surface
    else
        canopy_albedo[cell] = zero(T)
        pet_albedo[cell] = bare_surface
    end
end
