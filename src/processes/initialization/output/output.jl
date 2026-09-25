"""Time-series outputs for crop processes."""
mutable struct CropOutput{A, I, M}
    gpp::A                 # Daily gross primary production (gC m⁻² day⁻¹).
    npp::A                 # Daily net primary production after all plant respiration (gC m⁻² day⁻¹).
    lambda::A              # Optimal intercellular-to-ambient CO₂ ratio used by photosynthesis (0–1).
    potential_vcmax::A     # Water- and N-unlimited maximum carboxylation capacity (gC m⁻² day⁻¹).
    vcmax::A               # Realized maximum carboxylation capacity after limitations (gC m⁻² day⁻¹).
    nitrogen_limitation::A # Realized-to-potential `vcmax` ratio (0–1).
    respiration::A         # Total plant maintenance plus growth respiration (gC m⁻² day⁻¹).
    biomass::A             # Total live crop carbon stock (gC m⁻²).
    lai::A                 # Actual nonnegative leaf-area index (m² leaf m⁻² ground).
    storage_carbon::A      # Carbon in the harvestable storage organ (gC m⁻²).
    yield::A               # Annual harvested storage-organ carbon (gC m⁻² yr⁻¹).
    season_gpp::A          # Harvest-season cumulative gross primary production (gC m⁻²).
    season_lai_days::A     # Harvest-season cumulative LAI (m² leaf m⁻² ground day).
    season_length::A       # Active crop days from sowing through the day before harvest (day).
    season_water_deficit::A # Harvest-season cumulative crop water deficit (% day).
    season_evapotranspiration::A # Harvest-season total ET (mm).
    harvest_aboveground_carbon::A # Live above-ground crop carbon immediately before harvest (gC m⁻²).
    vegetation_carbon::M   # Daily leaf/root/pool/storage carbon stocks (gC m⁻²).
    vegetation_nitrogen::M # Daily leaf/root/pool/storage nitrogen contents (gN m⁻²).
    fphu::A                # Fraction of potential heat units accumulated (0–1+).
    water_deficit::A       # Daily crop water-deficit factor (0–100%).
    growing_mask::I        # Integer mask: one while a crop stand is active, otherwise zero.
end

"""Time-series outputs for soil stocks and fluxes."""
mutable struct SoilOutput{A, M}
    ecosystem_respiration::A   # Plant plus heterotrophic respiration (gC m⁻² day⁻¹).
    litter_carbon::M           # Carbon stocks in surface/incorporated/root litter (gC m⁻²).
    fast_carbon::M             # Fast soil-organic-carbon stock by layer (gC m⁻²).
    slow_carbon::M             # Slow soil-organic-carbon stock by layer (gC m⁻²).
    water_storage::M           # Liquid soil-water storage by layer (mm).
    litter_nitrogen::M         # Nitrogen stocks in the three litter classes (gN m⁻²).
    fast_nitrogen::M           # Fast soil-organic-nitrogen stock by layer (gN m⁻²).
    slow_nitrogen::M           # Slow soil-organic-nitrogen stock by layer (gN m⁻²).
    heterotrophic_respiration::A # Litter and soil respiration (gC m⁻² day⁻¹).
    evapotranspiration::A      # Total land-surface ET: canopy interception, litter/soil evaporation, and crop transpiration (mm day⁻¹).
end

"""Time-series outputs for climate forcing and potential evaporation."""
mutable struct ClimateOutput{A}
    equilibrium_evapotranspiration::A # Priestley–Taylor equilibrium evaporation demand (mm day⁻¹).
    precipitation::A                 # Daily precipitation forcing (mm day⁻¹).
    temperature::A                   # Daily near-surface air temperature forcing (°C).
end

"""Time-series crop calendar diagnostics."""
mutable struct CalendarOutput{I}
    harvesting_mask::I # Daily mask indicating the harvest window/condition.
    harvesting_year::I # Simulation year associated with each annual harvest record.
    harvest_date::I    # Day of year of the recorded annual harvest.
    sowing_event::I    # Daily one-day sowing event indicator (0/1).
    harvest_event::I   # Daily one-day harvest event indicator (0/1).
end

"""In-progress annual crop outputs retained until the calendar-year boundary."""
mutable struct AnnualOutputAccumulator{A, I}
    yield::A        # Harvested storage carbon accumulated in the current output year (gC m⁻²).
    harvest_date::I # Latest harvest day in the current output year (1–365; 0 if absent).
    season_gpp::A
    season_lai_days::A
    season_length::A
    season_water_deficit::A
    season_evapotranspiration::A
    harvest_aboveground_carbon::A
    active_gpp::A
    active_lai_days::A
    active_length::A
    active_water_deficit::A
    active_evapotranspiration::A
end

"""Process-grouped model output container."""
mutable struct Output{C, S, F, K, A}
    crop::C     # Crop daily and annual output time series.
    soil::S     # Soil stock and flux output time series.
    climate::F  # Selected climate-forcing output time series.
    calendar::K # Sowing and harvest calendar output time series.
    annual::A   # In-progress annual output records required before year-end emission.
end

const _DAILY_CROP_FLOAT_OUTPUT_FIELDS = (
    :gpp, :npp, :lambda, :potential_vcmax, :vcmax,
    :nitrogen_limitation, :respiration, :biomass, :lai,
    :storage_carbon, :fphu, :water_deficit,
)
const _DAILY_CROP_INTEGER_OUTPUT_FIELDS = (:growing_mask,)
const _DAILY_SOIL_FLOAT_OUTPUT_FIELDS = (
    :ecosystem_respiration, :heterotrophic_respiration, :evapotranspiration,
)
const _DAILY_CALENDAR_INTEGER_OUTPUT_FIELDS = (
    :harvesting_mask, :sowing_event, :harvest_event,
)
const _ANNUAL_CROP_FLOAT_OUTPUT_FIELDS = (
    :yield,
    :season_gpp,
    :season_lai_days,
    :season_length,
    :season_water_deficit,
    :season_evapotranspiration,
    :harvest_aboveground_carbon,
)
const _ANNUAL_CALENDAR_INTEGER_OUTPUT_FIELDS = (:harvest_date, :harvesting_year)

init_output(cell_size::Int, device; kwargs...) =
    init_output(Float32, cell_size, device; kwargs...)
function init_output(::Type{T},
                     cell_size::Int,
                     device;
                     vegc_pools::Int = 4,
                     litc_layers::Int = 3,
                     soil_layers::Int = 5) where {T <: AbstractFloat}
    # Output rows represent completed simulation steps only. Initial model
    # state lives in `crop`/`soil`; it is not a synthetic day-zero output.
    scalar_output() = device(zeros(T, 0, cell_size))
    integer_output() = device(zeros(Int32, 0, cell_size))

    crop = CropOutput(
        scalar_output(), scalar_output(), scalar_output(), scalar_output(),
        scalar_output(), scalar_output(), scalar_output(), scalar_output(),
        scalar_output(), scalar_output(), scalar_output(), scalar_output(),
        scalar_output(), scalar_output(), scalar_output(), scalar_output(),
        scalar_output(),
        device(zeros(T, 0, vegc_pools * cell_size)),
        device(zeros(T, 0, vegc_pools * cell_size)),
        scalar_output(), scalar_output(), integer_output(),
    )

    soil = SoilOutput(
        scalar_output(),
        device(zeros(T, 0, litc_layers * cell_size)),
        device(zeros(T, 0, soil_layers * cell_size)),
        device(zeros(T, 0, soil_layers * cell_size)),
        device(zeros(T, 0, soil_layers * cell_size)),
        device(zeros(T, 0, litc_layers * cell_size)),
        device(zeros(T, 0, soil_layers * cell_size)),
        device(zeros(T, 0, soil_layers * cell_size)),
        scalar_output(), scalar_output(),
    )

    climate = ClimateOutput(scalar_output(), scalar_output(), scalar_output())
    calendar = CalendarOutput(
        integer_output(), integer_output(), integer_output(),
        integer_output(), integer_output(),
    )
    annual = AnnualOutputAccumulator(
        device(zeros(T, cell_size)), device(zeros(Int32, cell_size)),
        device(zeros(T, cell_size)), device(zeros(T, cell_size)),
        device(zeros(T, cell_size)), device(zeros(T, cell_size)),
        device(zeros(T, cell_size)), device(zeros(T, cell_size)),
        device(zeros(T, cell_size)), device(zeros(T, cell_size)),
        device(zeros(T, cell_size)), device(zeros(T, cell_size)),
        device(zeros(T, cell_size)),
    )
    return Output(crop, soil, climate, calendar, annual)
end

"""Grow a backend array once for a simulation block, preserving existing rows."""
@kernel inbounds = true function copy_output_rows_kernel!(destination, source)
    row, column = @index(Global, NTuple)
    destination[row, column] = source[row, column]
end

@kernel inbounds = true function write_output_row_kernel!(destination, row, source)
    column = @index(Global)
    destination[row, column] = source[column]
end

function _extend_output_rows(array::AbstractMatrix, additional_rows::Integer)
    additional_rows <= 0 && return array
    old_rows, columns = size(array)
    extended = similar(array, old_rows + additional_rows, columns)
    fill!(extended, zero(eltype(extended)))
    old_rows > 0 && launch_custom!(
        copy_output_rows_kernel!, extended, (old_rows, columns), array,
    )
    return extended
end

function _reuse_output_rows(array::AbstractMatrix, rows::Integer)
    rows >= 0 || throw(ArgumentError("output rows must be non-negative"))
    if size(array, 1) == rows
        fill!(array, zero(eltype(array)))
        return array
    end
    resized = similar(array, rows, size(array, 2))
    fill!(resized, zero(eltype(resized)))
    return resized
end

"""
    prepare_output_block!(output, daily_rows, annual_rows)

Reserve all crop/calendar output rows once before a daily simulation block.
The returned indices point to the first newly allocated daily and annual rows.
"""
function prepare_output_block!(output::Output,
                               daily_rows::Integer,
                               annual_rows::Integer;
                               reuse::Bool = false,
                               selected::Union{Nothing, Set{Tuple{Symbol, Symbol}}} = nothing)
    !isnothing(selected) && !reuse && throw(ArgumentError(
        "selected output allocation requires reuse=true",
    ))
    first_daily_row = reuse ? 1 : size(output.crop.gpp, 1) + 1
    first_annual_row = reuse ? 1 : size(output.crop.yield, 1) + 1
    resize_rows = reuse ? _reuse_output_rows : _extend_output_rows

    for field in (_DAILY_CROP_FLOAT_OUTPUT_FIELDS..., _DAILY_CROP_INTEGER_OUTPUT_FIELDS...)
        rows = isnothing(selected) || (:crop, field) in selected ? daily_rows : 0
        setproperty!(
            output.crop,
            field,
            resize_rows(getproperty(output.crop, field), rows),
        )
    end
    for field in _DAILY_SOIL_FLOAT_OUTPUT_FIELDS
        rows = isnothing(selected) || (:soil, field) in selected ? daily_rows : 0
        setproperty!(
            output.soil,
            field,
            resize_rows(getproperty(output.soil, field), rows),
        )
    end
    for field in _DAILY_CALENDAR_INTEGER_OUTPUT_FIELDS
        rows = isnothing(selected) || (:calendar, field) in selected ? daily_rows : 0
        setproperty!(
            output.calendar,
            field,
            resize_rows(getproperty(output.calendar, field), rows),
        )
    end

    for field in _ANNUAL_CROP_FLOAT_OUTPUT_FIELDS
        rows = isnothing(selected) || (:crop, field) in selected ? annual_rows : 0
        setproperty!(output.crop, field, resize_rows(getproperty(output.crop, field), rows))
    end
    date_rows = isnothing(selected) || (:calendar, :harvest_date) in selected ? annual_rows : 0
    year_rows = isnothing(selected) || (:calendar, :harvesting_year) in selected ? annual_rows : 0
    output.calendar.harvest_date =
        resize_rows(output.calendar.harvest_date, date_rows)
    output.calendar.harvesting_year =
        resize_rows(output.calendar.harvesting_year, year_rows)

    return (; first_daily_row, first_annual_row)
end

@inline function _write_output_row!(destination::AbstractMatrix,
                                    row::Integer,
                                    source::AbstractVector)
    size(destination, 1) == 0 && return nothing
    launch_custom!(write_output_row_kernel!, destination, length(source), row, source)
    return nothing
end

function _append_output_row(array::AbstractMatrix, source::AbstractVector)
    row = size(array, 1) + 1
    extended = _extend_output_rows(array, 1)
    _write_output_row!(extended, row, source)
    return extended
end

"""
    record_ecosystem_flux_outputs!(output, crop, soil; output_row)

Write daily ecosystem flux diagnostics after the water and soil-carbon fluxes
are available. `ecosystem_respiration` is total ecosystem respiration and
`evapotranspiration` is total land-surface ET.
"""
function record_ecosystem_flux_outputs!(output::Output, crop, soil;
                                        output_row::Integer)
    plant_respiration = crop_fluxes(crop).carbon.respiration
    leaf_respiration = crop_fluxes(crop).carbon.leaf_respiration
    backend = KernelAbstractions.get_backend(plant_respiration)
    kernel = record_ecosystem_flux_outputs_kernel!(backend)
    kernel(
        output.soil.ecosystem_respiration,
        output.soil.heterotrophic_respiration,
        output.soil.evapotranspiration,
        plant_respiration,
        leaf_respiration,
        soil_carbon_fluxes(soil).heterotrophic_respiration,
        crop_fluxes(crop).water.interception,
        crop_fluxes(crop).water.transpiration_layer,
        soil_water_fluxes(soil).evaporation,
        soil_surface_litter_fluxes(soil).evaporation,
        output_row,
        size(crop_fluxes(crop).water.transpiration_layer, 1),
        ndrange = length(plant_respiration),
    )
    return nothing
end

"""Accumulate harvest-season diagnostics without feeding back into model state."""
function accumulate_season_process_diagnostics!(output::Output, crop, soil)
    water_flux = crop_fluxes(crop).water
    backend = KernelAbstractions.get_backend(crop_fluxes(crop).carbon.gross_assimilation)
    kernel = accumulate_season_process_diagnostics_kernel!(backend)
    kernel(
        output.annual.active_gpp,
        output.annual.active_lai_days,
        output.annual.active_length,
        output.annual.active_water_deficit,
        output.annual.active_evapotranspiration,
        crop_events(crop).sowing,
        crop_prognostic(crop).phenology.is_growing,
        crop_fluxes(crop).carbon.gross_assimilation,
        crop_canopy_auxiliary(crop).actual_lai,
        crop_stress_auxiliary(crop).water_deficit,
        water_flux.interception,
        water_flux.transpiration_layer,
        soil_water_fluxes(soil).evaporation,
        soil_surface_litter_fluxes(soil).evaporation,
        size(water_flux.transpiration_layer, 1),
        ndrange = length(crop_prognostic(crop).phenology.is_growing),
    )
    return nothing
end

@kernel inbounds = true function accumulate_season_process_diagnostics_kernel!(
    active_gpp::AbstractVector{T},
    active_lai_days::AbstractVector{T},
    active_length::AbstractVector{T},
    active_water_deficit::AbstractVector{T},
    active_evapotranspiration::AbstractVector{T},
    sowing_event::AbstractVector{S},
    is_growing::AbstractVector{S},
    gross_assimilation::AbstractVector{T},
    actual_lai::AbstractVector{T},
    water_deficit::AbstractVector{T},
    interception::AbstractVector{T},
    transpiration_layer::AbstractMatrix{T},
    soil_evaporation::AbstractMatrix{T},
    litter_evaporation::AbstractVector{T},
    layers::Integer,
) where {T <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    if sowing_event[cell] != 0
        active_gpp[cell] = zero(T)
        active_lai_days[cell] = zero(T)
        active_length[cell] = zero(T)
        active_water_deficit[cell] = zero(T)
        active_evapotranspiration[cell] = zero(T)
    end
    if is_growing[cell] != 0
        active_gpp[cell] += gross_assimilation[cell]
        active_lai_days[cell] += actual_lai[cell]
        active_length[cell] += one(T)
        active_water_deficit[cell] += water_deficit[cell]
        total_et = interception[cell] + litter_evaporation[cell]
        for layer in 1:layers
            total_et += transpiration_layer[layer, cell] + soil_evaporation[layer, cell]
        end
        active_evapotranspiration[cell] += total_et
    end
end

@kernel inbounds = true function record_ecosystem_flux_outputs_kernel!(
    ecosystem_respiration::AbstractMatrix{T},
    heterotrophic_respiration::AbstractMatrix{T},
    evapotranspiration::AbstractMatrix{T},
    plant_respiration::AbstractVector{T},
    leaf_respiration::AbstractVector{T},
    soil_heterotrophic_respiration::AbstractVector{T},
    canopy_interception::AbstractVector{T},
    transpiration_layer::AbstractMatrix{T},
    soil_evaporation::AbstractMatrix{T},
    litter_evaporation::AbstractVector{T},
    output_row::Integer,
    layers::Integer,
) where {T <: AbstractFloat}
    cell = @index(Global)
    heterotrophic = soil_heterotrophic_respiration[cell]
    if output_row <= size(ecosystem_respiration, 1)
        ecosystem_respiration[output_row, cell] =
            plant_respiration[cell] + leaf_respiration[cell] + heterotrophic
    end
    if output_row <= size(heterotrophic_respiration, 1)
        heterotrophic_respiration[output_row, cell] = heterotrophic
    end
    if output_row <= size(evapotranspiration, 1)
        total_et = canopy_interception[cell] + litter_evaporation[cell]
        for layer in 1:layers
            total_et += transpiration_layer[layer, cell] + soil_evaporation[layer, cell]
        end
        evapotranspiration[output_row, cell] = total_et
    end
end
