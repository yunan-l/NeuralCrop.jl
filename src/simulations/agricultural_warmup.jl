function _warmup_soil_snapshot(state::ModelState)
    carbon = soil_carbon_prognostic(state)
    nitrogen = soil_nitrogen_prognostic(state)
    water = soil_water_prognostic(state)
    host(values) = Array(values)

    litter_carbon = vec(sum(host(carbon.litter); dims = 1))
    litter_nitrogen = vec(sum(host(nitrogen.litter); dims = 1))
    fast_carbon = vec(sum(host(carbon.fast); dims = 1))
    slow_carbon = vec(sum(host(carbon.slow); dims = 1))
    fast_nitrogen = vec(sum(host(nitrogen.fast); dims = 1))
    slow_nitrogen = vec(sum(host(nitrogen.slow); dims = 1))
    mineral_nitrogen = vec(sum(
        host(nitrogen.nitrate) .+ host(nitrogen.ammonium); dims = 1,
    ))
    return (
        total_carbon = litter_carbon .+ fast_carbon .+ slow_carbon,
        total_nitrogen = litter_nitrogen .+ fast_nitrogen .+
            slow_nitrogen .+ mineral_nitrogen,
        litter_carbon = litter_carbon,
        litter_nitrogen = litter_nitrogen,
        fast_carbon = fast_carbon,
        fast_nitrogen = fast_nitrogen,
        slow_carbon = slow_carbon,
        slow_nitrogen = slow_nitrogen,
        mineral_nitrogen = mineral_nitrogen,
        soil_water = vec(sum(host(water.storage); dims = 1)),
    )
end

function _warmup_history_storage(initial, maximum_years::Integer)
    names = keys(initial)
    cells = length(first(initial))
    return NamedTuple{names}(map(
        name -> Matrix{eltype(getproperty(initial, name))}(undef, maximum_years, cells),
        names,
    ))
end

function _store_warmup_snapshot!(history, year::Integer, snapshot)
    for name in keys(history)
        getproperty(history, name)[year, :] .= getproperty(snapshot, name)
    end
    return history
end

_warmup_history_row(history, year::Integer) = map(
    values -> view(values, year, :), history,
)

_warmup_history_view(history, years::Integer) = map(
    values -> view(values, 1:years, :), history,
)

function _warmup_targets(state::ModelState)
    carbon = soil_carbon_prognostic(state)
    nitrogen = soil_nitrogen_prognostic(state)
    return (
        carbon = (fast = copy(carbon.fast), slow = copy(carbon.slow)),
        nitrogen = (
            fast = copy(nitrogen.fast),
            slow = copy(nitrogen.slow),
            nitrate = copy(nitrogen.nitrate),
            ammonium = copy(nitrogen.ammonium),
        ),
    )
end

function _warmup_c_shift_workspace(state::ModelState)
    decomposition = soil_decomposition_input(state)
    response = soil_decomposition_auxiliary(state).response
    response_sum = similar(response)
    fill!(response_sum, zero(eltype(response_sum)))
    return (
        response_sum,
        reference_fast = copy(decomposition.shift_fast),
        reference_slow = copy(decomposition.shift_slow),
    )
end

function _warmup_c_shift_report(workspace)
    fast = similar(workspace.reference_fast)
    slow = similar(workspace.reference_slow)
    equilibrated_c_shift!(
        fast, slow, workspace.response_sum, workspace.reference_fast,
        workspace.reference_slow,
    )
    return (
        fast = Array(fast),
        slow = Array(slow),
        response_sum = Array(workspace.response_sum),
    )
end

function _warmup_pool_allocation(state::ModelState, c_shift)
    carbon = Array(soil_carbon_prognostic(state).fast)
    carbon .+= Array(soil_carbon_prognostic(state).slow)
    nitrogen = Array(soil_nitrogen_prognostic(state).fast)
    nitrogen .+= Array(soil_nitrogen_prognostic(state).slow)
    fast_carbon = Array(soil_carbon_prognostic(state).fast)
    fast_nitrogen = Array(soil_nitrogen_prognostic(state).fast)
    return (
        fast_carbon_fraction = fast_carbon ./ max.(carbon, eps(eltype(carbon))),
        fast_nitrogen_fraction = fast_nitrogen ./ max.(nitrogen, eps(eltype(nitrogen))),
        c_shift_fast = c_shift.fast,
        c_shift_slow = c_shift.slow,
    )
end

@kernel inbounds = true function constrain_warmup_carbon_kernel!(
    fast::AbstractMatrix{T},
    slow::AbstractMatrix{T},
    target_fast::AbstractMatrix{T},
    target_slow::AbstractMatrix{T},
    correction::AbstractMatrix{T},
) where {T <: AbstractFloat}
    layer, cell = @index(Global, NTuple)
    target = target_fast[layer, cell] + target_slow[layer, cell]
    current = fast[layer, cell] + slow[layer, cell]
    correction[layer, cell] = target - current
    if current > eps(T)
        scale = target / current
        fast[layer, cell] *= scale
        slow[layer, cell] *= scale
    else
        fast[layer, cell] = target_fast[layer, cell]
        slow[layer, cell] = target_slow[layer, cell]
    end
end

@kernel inbounds = true function constrain_warmup_nitrogen_kernel!(
    fast::AbstractMatrix{T},
    slow::AbstractMatrix{T},
    nitrate::AbstractMatrix{T},
    ammonium::AbstractMatrix{T},
    target_fast::AbstractMatrix{T},
    target_slow::AbstractMatrix{T},
    target_nitrate::AbstractMatrix{T},
    target_ammonium::AbstractMatrix{T},
    correction::AbstractMatrix{T},
) where {T <: AbstractFloat}
    layer, cell = @index(Global, NTuple)
    target = target_fast[layer, cell] + target_slow[layer, cell] +
        target_nitrate[layer, cell] + target_ammonium[layer, cell]
    current = fast[layer, cell] + slow[layer, cell] +
        nitrate[layer, cell] + ammonium[layer, cell]
    correction[layer, cell] = target - current
    if current > eps(T)
        scale = target / current
        fast[layer, cell] *= scale
        slow[layer, cell] *= scale
        nitrate[layer, cell] *= scale
        ammonium[layer, cell] *= scale
    else
        fast[layer, cell] = target_fast[layer, cell]
        slow[layer, cell] = target_slow[layer, cell]
        nitrate[layer, cell] = target_nitrate[layer, cell]
        ammonium[layer, cell] = target_ammonium[layer, cell]
    end
end

function _apply_warmup_targets!(state::ModelState, targets)
    carbon = soil_carbon_prognostic(state)
    nitrogen = soil_nitrogen_prognostic(state)
    carbon_correction = similar(carbon.fast)
    nitrogen_correction = similar(nitrogen.fast)
    launch_2D!(
        constrain_warmup_carbon_kernel!, carbon.fast, carbon.slow,
        targets.carbon.fast, targets.carbon.slow, carbon_correction,
    )
    launch_2D!(
        constrain_warmup_nitrogen_kernel!, nitrogen.fast, nitrogen.slow,
        nitrogen.nitrate, nitrogen.ammonium, targets.nitrogen.fast,
        targets.nitrogen.slow, targets.nitrogen.nitrate,
        targets.nitrogen.ammonium, nitrogen_correction,
    )
    return (
        carbon = vec(sum(Array(carbon_correction); dims = 1)),
        nitrogen = vec(sum(Array(nitrogen_correction); dims = 1)),
    )
end

function _warmup_convergence_counts(
    consecutive,
    consecutive_years::Integer,
    convergence_reducer,
)
    converged_cells = count(>=(consecutive_years), consecutive)
    total_cells = length(consecutive)
    if convergence_reducer !== nothing
        converged_cells, total_cells = convergence_reducer(converged_cells, total_cells)
    end
    0 <= converged_cells <= total_cells || throw(ArgumentError(
        "warm-up convergence reduction returned invalid cell counts",
    ))
    total_cells > 0 || throw(ArgumentError(
        "warm-up convergence reduction returned no cells",
    ))
    return converged_cells, total_cells
end

function _warmup_convergence!(
    consecutive::AbstractVector{<:Integer},
    initial,
    previous,
    current,
    previous_correction,
    correction;
    relative_tolerance::Real,
    pool_fraction_tolerance::Real,
    consecutive_years::Integer,
    convergence_reducer = nothing,
)
    relative_change(now, before) = (now .- before) ./ max.(abs.(before), one(eltype(now)))
    pool_fraction(fast, slow) = fast ./ max.(fast .+ slow, eps(eltype(fast)))
    stable = abs.(relative_change(current.total_carbon, previous.total_carbon)) .<=
        relative_tolerance
    stable .&= abs.(relative_change(current.total_nitrogen, previous.total_nitrogen)) .<=
        relative_tolerance
    stable .&= abs.(
        pool_fraction(current.fast_carbon, current.slow_carbon) .-
        pool_fraction(previous.fast_carbon, previous.slow_carbon)
    ) .<= pool_fraction_tolerance
    stable .&= abs.(
        pool_fraction(current.fast_nitrogen, current.slow_nitrogen) .-
        pool_fraction(previous.fast_nitrogen, previous.slow_nitrogen)
    ) .<= pool_fraction_tolerance
    stable .&= abs.(correction.carbon .- previous_correction.carbon) ./
        max.(abs.(initial.total_carbon), one(eltype(correction.carbon))) .<= relative_tolerance
    stable .&= abs.(correction.nitrogen .- previous_correction.nitrogen) ./
        max.(abs.(initial.total_nitrogen), one(eltype(correction.nitrogen))) .<= relative_tolerance
    consecutive .= ifelse.(stable, consecutive .+ 1, 0)
    converged_cells, total_cells = _warmup_convergence_counts(
        consecutive, consecutive_years, convergence_reducer,
    )
    return converged_cells / total_cells
end

function _warmup_options(
    years,
    maximum_years,
    consecutive_years,
    relative_tolerance,
    pool_fraction_tolerance,
    required_converged_fraction,
)
    years > 0 || throw(ArgumentError("warm-up years must be positive"))
    maximum_years >= years || throw(ArgumentError(
        "maximum warm-up years must be at least the minimum years",
    ))
    consecutive_years > 0 || throw(ArgumentError("consecutive warm-up years must be positive"))
    relative_tolerance > 0 || throw(ArgumentError("warm-up relative tolerance must be positive"))
    pool_fraction_tolerance > 0 || throw(ArgumentError(
        "warm-up pool-fraction tolerance must be positive",
    ))
    0 < required_converged_fraction <= 1 || throw(ArgumentError(
        "required converged fraction must be in (0, 1]",
    ))
    return nothing
end

"""Summarize initial, transient, and late-year C/N drift from a warm-up report."""
function agricultural_warmup_drift(
    report;
    initial_fast_fraction_threshold::Real = 0.10,
    late_fast_fraction_threshold::Real = 0.01,
    late_total_threshold::Real = 0.01,
)
    report.years >= 2 || throw(ArgumentError("C/N drift requires at least two warm-up years"))
    series(name) = vcat(
        sum(getproperty(report.initial_soil, name)),
        vec(sum(getproperty(report.soil, name); dims = 2)),
    )
    total_carbon = series(:total_carbon)
    total_nitrogen = series(:total_nitrogen)
    fast_carbon = series(:fast_carbon)
    slow_carbon = series(:slow_carbon)
    fast_fraction = fast_carbon ./ (fast_carbon .+ slow_carbon)
    initial_fast_shift = fast_fraction[2] - fast_fraction[1]
    forcing_years = hasproperty(report, :forcing_years) ? report.forcing_years : 1
    comparison_index = max(1, length(total_carbon) - forcing_years)
    late_fast_shift = fast_fraction[end] - fast_fraction[comparison_index]
    late_carbon = (total_carbon[end] - total_carbon[comparison_index]) /
        max(abs(total_carbon[comparison_index]), eps(eltype(total_carbon)))
    late_nitrogen = (total_nitrogen[end] - total_nitrogen[comparison_index]) /
        max(abs(total_nitrogen[comparison_index]), eps(eltype(total_nitrogen)))
    relative_values(current, previous) = (current .- previous) ./
        max.(abs.(previous), eps(eltype(previous)))
    distribution(values) = (
        minimum = minimum(values),
        p05 = quantile(values, 0.05),
        p25 = quantile(values, 0.25),
        median = median(values),
        p75 = quantile(values, 0.75),
        p95 = quantile(values, 0.95),
        maximum = maximum(values),
    )
    initial_carbon = report.initial_soil.total_carbon
    initial_nitrogen = report.initial_soil.total_nitrogen
    final_carbon = vec(report.soil.total_carbon[end, :])
    final_nitrogen = vec(report.soil.total_nitrogen[end, :])
    comparison_year = max(0, report.years - forcing_years)
    previous_carbon = comparison_year == 0 ? initial_carbon :
        vec(report.soil.total_carbon[comparison_year, :])
    previous_nitrogen = comparison_year == 0 ? initial_nitrogen :
        vec(report.soil.total_nitrogen[comparison_year, :])
    initial_cell_fast_fraction = report.initial_soil.fast_carbon ./
        (report.initial_soil.fast_carbon .+ report.initial_soil.slow_carbon)
    final_cell_fast_fraction = vec(report.soil.fast_carbon[end, :]) ./
        (vec(report.soil.fast_carbon[end, :]) .+ vec(report.soil.slow_carbon[end, :]))
    previous_fast_carbon = comparison_year == 0 ? report.initial_soil.fast_carbon :
        vec(report.soil.fast_carbon[comparison_year, :])
    previous_slow_carbon = comparison_year == 0 ? report.initial_soil.slow_carbon :
        vec(report.soil.slow_carbon[comparison_year, :])
    previous_cell_fast_fraction = previous_fast_carbon ./
        (previous_fast_carbon .+ previous_slow_carbon)
    cell_initial_carbon = relative_values(final_carbon, initial_carbon)
    cell_initial_nitrogen = relative_values(final_nitrogen, initial_nitrogen)
    cell_late_carbon = relative_values(final_carbon, previous_carbon)
    cell_late_nitrogen = relative_values(final_nitrogen, previous_nitrogen)
    cell_initial_fast_shift = final_cell_fast_fraction .- initial_cell_fast_fraction
    cell_late_fast_shift = final_cell_fast_fraction .- previous_cell_fast_fraction
    carbon_correction = hasproperty(report, :target_correction) ?
        report.target_correction.carbon : zeros(eltype(final_carbon), report.years, length(final_carbon))
    nitrogen_correction = hasproperty(report, :target_correction) ?
        report.target_correction.nitrogen : zeros(eltype(final_nitrogen), report.years, length(final_nitrogen))
    final_carbon_correction = vec(carbon_correction[end, :])
    final_nitrogen_correction = vec(nitrogen_correction[end, :])
    relative_carbon_correction = final_carbon_correction ./
        max.(abs.(initial_carbon), one(eltype(final_carbon_correction)))
    relative_nitrogen_correction = final_nitrogen_correction ./
        max.(abs.(initial_nitrogen), one(eltype(final_nitrogen_correction)))
    cell_review = (abs.(cell_initial_fast_shift) .> initial_fast_fraction_threshold) .|
        (abs.(cell_late_fast_shift) .> late_fast_fraction_threshold) .|
        (abs.(cell_late_carbon) .> late_total_threshold) .|
        (abs.(cell_late_nitrogen) .> late_total_threshold)
    review = abs(initial_fast_shift) > initial_fast_fraction_threshold ||
        abs(late_fast_shift) > late_fast_fraction_threshold ||
        abs(late_carbon) > late_total_threshold ||
        abs(late_nitrogen) > late_total_threshold
    target_constrained = hasproperty(report, :target_constrained) && report.target_constrained
    recommendation = if target_constrained
        report.converged ? :target_constrained_converged : :target_constrained_maximum_years
    else
        review ? :review_pool_allocation : :retain_40_60
    end
    return (
        total_carbon,
        total_nitrogen,
        fast_carbon_fraction = fast_fraction,
        initial_fast_fraction_shift = initial_fast_shift,
        late_fast_fraction_shift = late_fast_shift,
        late_carbon_relative_change = late_carbon,
        late_nitrogen_relative_change = late_nitrogen,
        convergence = (
            target_constrained,
            comparison_lag_years = forcing_years,
            converged = hasproperty(report, :converged) ? report.converged : false,
            actual_years = report.years,
            converged_cell_fraction = hasproperty(report, :converged_cell_fraction) ?
                report.converged_cell_fraction : 0.0,
            unconverged_cells = hasproperty(report, :unconverged_cells) ?
                report.unconverged_cells : length(cell_review),
        ),
        target_correction = (
            annual_carbon = vec(sum(carbon_correction; dims = 2)),
            annual_nitrogen = vec(sum(nitrogen_correction; dims = 2)),
            final_carbon_relative = distribution(relative_carbon_correction),
            final_nitrogen_relative = distribution(relative_nitrogen_correction),
        ),
        spatial = (
            initial_to_final_carbon = distribution(cell_initial_carbon),
            initial_to_final_nitrogen = distribution(cell_initial_nitrogen),
            same_phase_to_final_carbon = distribution(cell_late_carbon),
            same_phase_to_final_nitrogen = distribution(cell_late_nitrogen),
            initial_to_final_fast_fraction = distribution(cell_initial_fast_shift),
            same_phase_to_final_fast_fraction = distribution(cell_late_fast_shift),
            review_cell_fraction = count(cell_review) / length(cell_review),
            cell_count = length(cell_review),
        ),
        recommendation,
    )
end

"""
    agricultural_warmup!(simulation, climate; years=10, maximum_years=years,
                         target_constrained=false)

Cycle complete agricultural forcing years before the production run. The
forcing must contain one or more whole 365-day years; multi-year blocks cycle
in order when `years` exceeds the available history. Crop, soil, water,
thermal, carbon, and nitrogen processes are active, but
production outputs, balance ledgers, and `simulation.simulated_days` are left
untouched. The warmed prognostic state is retained.

When `target_constrained=true`, the initial mineral-soil C and total-N stocks
are restored after each year while litter remains unconstrained. Convergence
compares states at the same phase of the forcing cycle. After the minimum
`years`, annual cycling continues until the requested fraction of cells has met
the C/N, pool-fraction, and target-correction tolerances for
`consecutive_years`, or `maximum_years` is reached.

The returned report contains one host-side row per completed warm-up year and
one column per active cell.
"""
function agricultural_warmup!(
    simulation::CropSimulation,
    climate::NamedTuple;
    years::Integer = 10,
    maximum_years::Integer = years,
    target_constrained::Bool = false,
    consecutive_years::Integer = 3,
    relative_tolerance::Real = 0.01,
    pool_fraction_tolerance::Real = 0.01,
    required_converged_fraction::Real = 1.0,
    convergence_reducer = nothing,
)
    _warmup_options(
        years, maximum_years, consecutive_years, relative_tolerance,
        pool_fraction_tolerance, required_converged_fraction,
    )
    simulation.simulated_days == 0 || throw(ArgumentError(
        "agricultural warm-up must run before the production simulation",
    ))
    _output_timeseries_empty(simulation.output) || throw(ArgumentError(
        "agricultural warm-up requires empty production output",
    ))

    prepared_climate = _prepare_climate(simulation, climate)
    climate_days = size(prepared_climate.temp, 1)
    climate_days > 0 && climate_days % 365 == 0 || throw(ArgumentError(
        "agricultural warm-up requires complete 365-day climate years, got $climate_days rows",
    ))
    cells = length(simulation.config.execution.domain.indices)
    size(prepared_climate.temp, 2) == cells || throw(DimensionMismatch(
        "warm-up climate has $(size(prepared_climate.temp, 2)) cells; simulation has $cells",
    ))

    initial_soil = _warmup_soil_snapshot(simulation.state)
    history = _warmup_history_storage(initial_soil, maximum_years)
    correction_type = eltype(initial_soil.total_carbon)
    carbon_correction = zeros(correction_type, maximum_years, cells)
    nitrogen_correction = zeros(correction_type, maximum_years, cells)
    converged_fraction = zeros(Float64, maximum_years)
    targets = target_constrained ? _warmup_targets(simulation.state) : nothing
    c_shift_workspace = _warmup_c_shift_workspace(simulation.state)
    consecutive = zeros(Int, cells)
    actual_years = 0
    converged = false
    no_output = Set{Tuple{Symbol, Symbol}}()
    common = (
        irrigation = simulation.config.irrigation,
        manure = simulation.config.manure,
        fertilizer = simulation.config.fertilizer,
        with_tillage = simulation.config.with_tillage,
        crop_resp_fix = simulation.config.crop_resp_fix,
        nitrogen_limit_vcmax = simulation.config.nitrogen_limit_vcmax,
        reuse_output = true,
        selected_output = no_output,
        c_shift_response_sum = c_shift_workspace.response_sum,
    )

    forcing_years = div(climate_days, 365)
    for year in 1:maximum_years
        forcing_year = mod(year - 1, forcing_years)
        start_day = 365 * forcing_year + 1
        end_day = start_day + 364
        _daily_crop!(
            pathway_value(simulation.processes.pathway),
            start_day, end_day, simulation.processes, prepared_climate, simulation.state;
            simulation_day_offset = 365 * (year - 1) + 1 - start_day,
            # LPJmL updates V_req while crop dates/PHU are being established,
            # then freezes it for prescribed fixed-date production. Day one of
            # forcing year N+1 closes year N, hence the extra year here.
            update_vernalization_requirement =
                !simulation.config.freeze_vernalization_requirement || year <= forcing_years + 1,
            common...,
        )
        correction = if target_constrained
            _apply_warmup_targets!(simulation.state, targets)
        else
            (carbon = zeros(correction_type, cells), nitrogen = zeros(correction_type, cells))
        end
        carbon_correction[year, :] .= correction.carbon
        nitrogen_correction[year, :] .= correction.nitrogen
        current_soil = _warmup_soil_snapshot(simulation.state)
        _store_warmup_snapshot!(history, year, current_soil)
        # A persistent external input such as atmospheric N deposition can
        # require a stable, non-zero annual target correction. Constrained
        # convergence therefore compares corrections at the same forcing
        # phase instead of requiring the correction itself to vanish.
        if year >= forcing_years && (!target_constrained || year > forcing_years)
            previous_soil = year == forcing_years ? initial_soil :
                _warmup_history_row(history, year - forcing_years)
            previous_correction = target_constrained ? (
                carbon = view(carbon_correction, year - forcing_years, :),
                nitrogen = view(nitrogen_correction, year - forcing_years, :),
            ) : correction
            converged_fraction[year] = _warmup_convergence!(
                consecutive, initial_soil, previous_soil, current_soil,
                previous_correction, correction;
                relative_tolerance, pool_fraction_tolerance, consecutive_years,
                convergence_reducer,
            )
        end
        actual_years = year
        if year >= years && converged_fraction[year] >= required_converged_fraction
            converged = true
            break
        end
    end

    # `update_climbuf!` normally closes a year on the following day 1. The
    # production clock intentionally restarts at day 1, so close the final
    # warm-up year explicitly before returning.
    annual_climbuf!(
        simulation.climbuf.atemp, simulation.climbuf, simulation.cft;
        update_vernalization_requirement = !simulation.config.freeze_vernalization_requirement,
    )
    clear_output_timeseries!(simulation.output)
    calibrated_c_shift = _warmup_c_shift_report(c_shift_workspace)
    converged_cells, total_cells = _warmup_convergence_counts(
        consecutive, consecutive_years, convergence_reducer,
    )
    return (
        years = actual_years,
        days = 365 * actual_years,
        forcing_years = forcing_years,
        minimum_years = Int(years),
        maximum_years = Int(maximum_years),
        target_constrained,
        consecutive_years = Int(consecutive_years),
        relative_tolerance = Float64(relative_tolerance),
        pool_fraction_tolerance = Float64(pool_fraction_tolerance),
        required_converged_fraction = Float64(required_converged_fraction),
        converged,
        converged_cell_fraction = converged_fraction[actual_years],
        unconverged_cells = total_cells - converged_cells,
        consecutive_stable_years = consecutive,
        initial_soil = initial_soil,
        soil = _warmup_history_view(history, actual_years),
        target_correction = (
            carbon = carbon_correction[1:actual_years, :],
            nitrogen = nitrogen_correction[1:actual_years, :],
        ),
        converged_fraction = converged_fraction[1:actual_years],
        calibrated_c_shift,
        calibrated_pool_allocation = _warmup_pool_allocation(
            simulation.state, calibrated_c_shift,
        ),
    )
end

"""
    agricultural_warmup!(simulation, climate_blocks; years=10,
                         maximum_years=years, target_constrained=false)

Run the agricultural warm-up from restartable climate blocks without joining
the forcing into one in-memory array. Block boundaries may occur anywhere
within a 365-day forcing year.
"""
function agricultural_warmup!(
    simulation::CropSimulation,
    climate_blocks::AbstractVector;
    years::Integer = 10,
    maximum_years::Integer = years,
    target_constrained::Bool = false,
    consecutive_years::Integer = 3,
    relative_tolerance::Real = 0.01,
    pool_fraction_tolerance::Real = 0.01,
    required_converged_fraction::Real = 1.0,
    management_blocks::Union{Nothing, AbstractVector} = nothing,
    convergence_reducer = nothing,
)
    _warmup_options(
        years, maximum_years, consecutive_years, relative_tolerance,
        pool_fraction_tolerance, required_converged_fraction,
    )
    simulation.simulated_days == 0 || throw(ArgumentError(
        "agricultural warm-up must run before the production simulation",
    ))
    _output_timeseries_empty(simulation.output) || throw(ArgumentError(
        "agricultural warm-up requires empty production output",
    ))
    isempty(climate_blocks) && throw(ArgumentError(
        "agricultural warm-up requires at least one climate block",
    ))

    block_days = map(eachindex(climate_blocks)) do index
        block = climate_blocks[index]
        hasproperty(block, :temp) || throw(ArgumentError(
            "warm-up climate block $index has no temp field",
        ))
        size(block.temp, 1)
    end
    climate_days = sum(block_days)
    climate_days > 0 && climate_days % 365 == 0 || throw(ArgumentError(
        "agricultural warm-up requires complete 365-day climate years, got $climate_days rows",
    ))

    cells = length(simulation.config.execution.domain.indices)
    initial_soil = _warmup_soil_snapshot(simulation.state)
    history = _warmup_history_storage(initial_soil, maximum_years)
    correction_type = eltype(initial_soil.total_carbon)
    carbon_correction = zeros(correction_type, maximum_years, cells)
    nitrogen_correction = zeros(correction_type, maximum_years, cells)
    converged_fraction = zeros(Float64, maximum_years)
    targets = target_constrained ? _warmup_targets(simulation.state) : nothing
    c_shift_workspace = _warmup_c_shift_workspace(simulation.state)
    consecutive = zeros(Int, cells)
    actual_years = 0
    converged = false
    no_output = Set{Tuple{Symbol, Symbol}}()
    common = (
        irrigation = simulation.config.irrigation,
        manure = simulation.config.manure,
        fertilizer = simulation.config.fertilizer,
        with_tillage = simulation.config.with_tillage,
        crop_resp_fix = simulation.config.crop_resp_fix,
        nitrogen_limit_vcmax = simulation.config.nitrogen_limit_vcmax,
        reuse_output = true,
        selected_output = no_output,
        c_shift_response_sum = c_shift_workspace.response_sum,
    )

    forcing_years = div(climate_days, 365)
    management_blocks === nothing || length(management_blocks) == forcing_years ||
        throw(DimensionMismatch(
            "management_blocks must contain one row for each of the $forcing_years forcing years",
        ))
    block_ends = cumsum(block_days)
    warmup_day = 0
    for year in 1:maximum_years
        forcing_year = mod(year - 1, forcing_years)
        forcing_start = 365 * forcing_year + 1
        forcing_end = forcing_start + 364
        prescribed = management_blocks === nothing ?
            (prescribed_phu = nothing, prescribed_winter_type = nothing) :
            _prepare_annual_management!(
                simulation, management_blocks[forcing_year + 1],
            )

        for index in eachindex(climate_blocks)
            block_start = block_ends[index] - block_days[index] + 1
            block_end = block_ends[index]
            segment_start = max(forcing_start, block_start)
            segment_end = min(forcing_end, block_end)
            segment_start <= segment_end || continue

            prepared_climate = _prepare_climate(simulation, climate_blocks[index])
            size(prepared_climate.temp, 2) == cells || throw(DimensionMismatch(
                "warm-up climate block $index has $(size(prepared_climate.temp, 2)) cells; simulation has $cells",
            ))
            local_start = segment_start - block_start + 1
            local_end = segment_end - block_start + 1
            _daily_crop!(
                pathway_value(simulation.processes.pathway),
                local_start, local_end, simulation.processes, prepared_climate,
                simulation.state;
                simulation_day_offset = warmup_day + 1 - local_start,
                update_vernalization_requirement =
                    !simulation.config.freeze_vernalization_requirement || year <= forcing_years + 1,
                prescribed...,
                common...,
            )
            warmup_day += local_end - local_start + 1
        end
        correction = if target_constrained
            _apply_warmup_targets!(simulation.state, targets)
        else
            (carbon = zeros(correction_type, cells), nitrogen = zeros(correction_type, cells))
        end
        carbon_correction[year, :] .= correction.carbon
        nitrogen_correction[year, :] .= correction.nitrogen
        current_soil = _warmup_soil_snapshot(simulation.state)
        _store_warmup_snapshot!(history, year, current_soil)
        if year >= forcing_years && (!target_constrained || year > forcing_years)
            previous_soil = year == forcing_years ? initial_soil :
                _warmup_history_row(history, year - forcing_years)
            previous_correction = target_constrained ? (
                carbon = view(carbon_correction, year - forcing_years, :),
                nitrogen = view(nitrogen_correction, year - forcing_years, :),
            ) : correction
            converged_fraction[year] = _warmup_convergence!(
                consecutive, initial_soil, previous_soil, current_soil,
                previous_correction, correction;
                relative_tolerance, pool_fraction_tolerance, consecutive_years,
                convergence_reducer,
            )
        end
        actual_years = year
        if year >= years && converged_fraction[year] >= required_converged_fraction
            converged = true
            break
        end
    end

    annual_climbuf!(
        simulation.climbuf.atemp, simulation.climbuf, simulation.cft;
        update_vernalization_requirement = !simulation.config.freeze_vernalization_requirement,
    )
    clear_output_timeseries!(simulation.output)
    calibrated_c_shift = _warmup_c_shift_report(c_shift_workspace)
    converged_cells, total_cells = _warmup_convergence_counts(
        consecutive, consecutive_years, convergence_reducer,
    )
    return (
        years = actual_years,
        days = 365 * actual_years,
        forcing_years = forcing_years,
        minimum_years = Int(years),
        maximum_years = Int(maximum_years),
        target_constrained,
        consecutive_years = Int(consecutive_years),
        relative_tolerance = Float64(relative_tolerance),
        pool_fraction_tolerance = Float64(pool_fraction_tolerance),
        required_converged_fraction = Float64(required_converged_fraction),
        converged,
        converged_cell_fraction = converged_fraction[actual_years],
        unconverged_cells = total_cells - converged_cells,
        consecutive_stable_years = consecutive,
        initial_soil = initial_soil,
        soil = _warmup_history_view(history, actual_years),
        target_correction = (
            carbon = carbon_correction[1:actual_years, :],
            nitrogen = nitrogen_correction[1:actual_years, :],
        ),
        converged_fraction = converged_fraction[1:actual_years],
        calibrated_c_shift,
        calibrated_pool_allocation = _warmup_pool_allocation(
            simulation.state, calibrated_c_shift,
        ),
    )
end
