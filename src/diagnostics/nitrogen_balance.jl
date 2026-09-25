"""
    NitrogenBalance

Optional daily nitrogen-budget diagnostics. All fields have shape
`(number_of_days, number_of_cells)` and use the model nitrogen unit (g N m⁻²).

The tracked system stock is plant total N plus soil mineral and organic N.
Organ N pools are derived partitions of plant total N and are not counted
again. Positive `residual` denotes nitrogen entering the tracked system
without appearing in final storage or a recorded boundary loss.
"""
mutable struct NitrogenBalance{M <: AbstractArray{<:AbstractFloat}}
    plant_before::M                # Plant N stock at start of day (gN m⁻²).
    plant_after::M                 # Plant N stock at end of day (gN m⁻²).
    mineral_before::M              # Soil NO₃ plus NH₄ stock at start of day (gN m⁻²).
    mineral_after::M               # Soil NO₃ plus NH₄ stock at end of day (gN m⁻²).
    organic_before::M              # Litter plus organic-soil N at start of day (gN m⁻²).
    organic_after::M               # Litter plus organic-soil N at end of day (gN m⁻²).
    total_before::M                # Tracked ecosystem N stock at start of day (gN m⁻²).
    total_after::M                 # Tracked ecosystem N stock at end of day (gN m⁻²).
    root_uptake::M                 # Soil mineral N transferred to crop roots (gN m⁻² day⁻¹).
    seed_input::M                  # Seed nitrogen boundary input (gN m⁻² day⁻¹).
    prescribed_fertilizer_input::M # Scheduled mineral-fertilizer input (gN m⁻² day⁻¹).
    prescribed_manure_input::M     # Scheduled manure-N input (gN m⁻² day⁻¹).
    automatic_fertilizer_input::M  # Demand-driven mineral-N input (gN m⁻² day⁻¹).
    atmospheric_deposition::M      # Nitrate plus ammonium deposition input (gN m⁻² day⁻¹).
    harvest_export::M              # Harvested nitrogen leaving the system (gN m⁻² day⁻¹).
    mineralization::M              # Organic N mineralized today (gN m⁻² day⁻¹).
    immobilization::M              # Mineral N immobilized today (gN m⁻² day⁻¹).
    nitrification::M               # NH₄ consumed by nitrification (gN m⁻² day⁻¹).
    n2o_nitrification::M           # N₂O-N loss from nitrification (gN m⁻² day⁻¹).
    denitrification::M             # NO₃ consumed by denitrification (gN m⁻² day⁻¹).
    n2o_denitrification::M         # N₂O-N loss from denitrification (gN m⁻² day⁻¹).
    n2_denitrification::M          # N₂-N loss from denitrification (gN m⁻² day⁻¹).
    volatilization::M              # NH₃-N volatilization loss (gN m⁻² day⁻¹).
    gaseous_loss::M                # Sum of all recorded gaseous N losses (gN m⁻² day⁻¹).
    leaching_loss::M               # Mineral-N leaching loss (gN m⁻² day⁻¹).
    residual::M                    # Absolute daily nitrogen-budget closure error (gN m⁻²).
    relative_residual::M           # Nitrogen residual normalized by daily budget magnitude.
end

function init_nitrogen_balance(number_of_days::Integer,
                               number_of_cells::Integer,
                               device = identity;
                               T::Type{<:AbstractFloat} = Float32)
    allocate() = device(zeros(T, number_of_days, number_of_cells))
    return NitrogenBalance(
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(),
    )
end

function record_nitrogen_balance_start!(balance::NitrogenBalance,
                                        day_index::Integer,
                                        crop,
                                        soil)
    @views begin
        balance.plant_before[day_index, :] .= crop_prognostic(crop).nitrogen.total
        balance.mineral_before[day_index, :] .= vec(sum(
            soil_nitrogen_prognostic(soil).nitrate .+ soil_nitrogen_prognostic(soil).ammonium; dims = 1,
        ))
        balance.organic_before[day_index, :] .= vec(sum(
            soil_nitrogen_prognostic(soil).litter; dims = 1,
        )) .+ vec(sum(soil_nitrogen_prognostic(soil).fast .+ soil_nitrogen_prognostic(soil).slow; dims = 1))
        balance.total_before[day_index, :] .=
            balance.plant_before[day_index, :] .+
            balance.mineral_before[day_index, :] .+
            balance.organic_before[day_index, :]
    end
    return nothing
end

function record_nitrogen_balance_end!(balance::NitrogenBalance,
                                      day_index::Integer,
                                      crop,
                                      soil;
                                      nitrate_deposition = nothing,
                                      ammonium_deposition = nothing)
    @views begin
        balance.plant_after[day_index, :] .= crop_prognostic(crop).nitrogen.total
        balance.mineral_after[day_index, :] .= vec(sum(
            soil_nitrogen_prognostic(soil).nitrate .+ soil_nitrogen_prognostic(soil).ammonium; dims = 1,
        ))
        balance.organic_after[day_index, :] .= vec(sum(
            soil_nitrogen_prognostic(soil).litter; dims = 1,
        )) .+ vec(sum(soil_nitrogen_prognostic(soil).fast .+ soil_nitrogen_prognostic(soil).slow; dims = 1))
        balance.total_after[day_index, :] .=
            balance.plant_after[day_index, :] .+
            balance.mineral_after[day_index, :] .+
            balance.organic_after[day_index, :]

        balance.root_uptake[day_index, :] .=
            crop_fluxes(crop).nitrogen.uptake .- crop_fluxes(crop).nitrogen.auto_fertilizer
        balance.seed_input[day_index, :] .= crop_fluxes(crop).nitrogen.seed_input
        balance.prescribed_fertilizer_input[day_index, :] .=
            crop_fluxes(crop).nitrogen.prescribed_fertilizer_input
        balance.prescribed_manure_input[day_index, :] .=
            crop_fluxes(crop).nitrogen.prescribed_manure_input
        balance.automatic_fertilizer_input[day_index, :] .=
            crop_fluxes(crop).nitrogen.auto_fertilizer
        if isnothing(nitrate_deposition)
            balance.atmospheric_deposition[day_index, :] .= 0
        else
            balance.atmospheric_deposition[day_index, :] .=
                nitrate_deposition .+ ammonium_deposition
        end
        balance.harvest_export[day_index, :] .= crop_fluxes(crop).nitrogen.harvest_export
        balance.mineralization[day_index, :] .=
            vec(sum(soil_nitrogen_fluxes(soil).mineralization; dims = 1))
        balance.immobilization[day_index, :] .=
            vec(sum(soil_nitrogen_fluxes(soil).immobilization; dims = 1))
        balance.nitrification[day_index, :] .=
            vec(sum(soil_nitrogen_fluxes(soil).nitrification; dims = 1))
        balance.n2o_nitrification[day_index, :] .=
            vec(sum(soil_nitrogen_fluxes(soil).n2o_nitrification; dims = 1))
        balance.denitrification[day_index, :] .=
            vec(sum(soil_nitrogen_fluxes(soil).denitrification; dims = 1))
        balance.n2o_denitrification[day_index, :] .=
            vec(sum(soil_nitrogen_fluxes(soil).n2o_denitrification; dims = 1))
        balance.n2_denitrification[day_index, :] .=
            vec(sum(soil_nitrogen_fluxes(soil).n2_denitrification; dims = 1))
        balance.volatilization[day_index, :] .= soil_nitrogen_fluxes(soil).volatilization
        balance.gaseous_loss[day_index, :] .=
            balance.n2o_nitrification[day_index, :] .+
            balance.n2o_denitrification[day_index, :] .+
            balance.n2_denitrification[day_index, :] .+
            balance.volatilization[day_index, :]
        balance.leaching_loss[day_index, :] .= soil_nitrogen_fluxes(soil).leaching

        balance.residual[day_index, :] .=
            balance.total_before[day_index, :] .+
            balance.seed_input[day_index, :] .+
            balance.prescribed_fertilizer_input[day_index, :] .+
            balance.prescribed_manure_input[day_index, :] .+
            balance.automatic_fertilizer_input[day_index, :] .+
            balance.atmospheric_deposition[day_index, :] .-
            balance.harvest_export[day_index, :] .-
            balance.gaseous_loss[day_index, :] .-
            balance.leaching_loss[day_index, :] .-
            balance.total_after[day_index, :]
        balance.relative_residual[day_index, :] .=
            balance.residual[day_index, :] ./ max.(
                abs.(balance.total_before[day_index, :]) .+
                balance.seed_input[day_index, :] .+
                balance.prescribed_fertilizer_input[day_index, :] .+
                balance.prescribed_manure_input[day_index, :] .+
                balance.automatic_fertilizer_input[day_index, :] .+
                balance.atmospheric_deposition[day_index, :],
                eps(eltype(balance.residual)),
            )
    end
    return nothing
end
