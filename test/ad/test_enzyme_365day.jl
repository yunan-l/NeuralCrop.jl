function _fixed_365_context(template, days)
    values = _seasonal_production_values(template, days)
    context = ADSeasonContext(
        trues(length(days)),
        (gpp = trues(length(days)), reco = trues(length(days)), et = trues(length(days))),
        (
            gpp = values.gpp .+ 0.1f0,
            reco = values.reco .+ 0.05f0,
            et = values.et .+ 0.05f0,
        ),
        (gpp = 1.0f0, reco = 1.0f0, et = 1.0f0),
    )
    return values, context
end

@testset "Enzyme fixed-shape 365-day season" begin
    # Establish the crop and fixed discrete event state first; the AD window
    # then covers exactly 365 continuous days after that boundary.
    template = _daily_transition_fixture(375, 10)
    days = Tuple(11:375)
    production_values, context = _fixed_365_context(template, days)
    parameter_names = (:gmin, :b)
    theta = Float32[cft1.gmin, cft1.b]
    direction = Float32[1, -1] ./ sqrt(2.0f0)

    adapter_value = _seasonal_loss_value(template, theta, parameter_names, days, context)
    production_value = _expected_seasonal_loss(production_values, context)
    @test isfinite(adapter_value)
    @test isfinite(production_value)
    # The production daily driver also executes cultivate/harvest and other
    # discrete events. This fixed-event AD fixture intentionally excludes
    # those events, so production loss is diagnostic rather than a parity
    # assertion for the continuous trajectory.
    @info "365-day fixed-event loss" ad_loss = adapter_value production_loss = production_value

    state_factory() = begin
        state = deepcopy(template.state)
        return state, enzyme_zero_tangent(state)
    end
    ad = enzyme_forward_directional(
        enzyme_seasonal_loss,
        theta,
        direction,
        state_factory,
        cft1,
        template.global_parameters,
        template.climate,
        parameter_names,
        days,
        template.layer_depth,
        context;
        return_primal = true,
    )
    @test isfinite(ad.primal)
    @test isfinite(ad.directional)
    @test ad.primal ≈ adapter_value rtol = 1.0f-6 atol = 1.0f-7

    step = 1.0f-3
    finite_difference = (
        _seasonal_loss_value(
            template, theta .+ step .* direction, parameter_names, days, context,
        ) -
        _seasonal_loss_value(
            template, theta .- step .* direction, parameter_names, days, context,
        )
    ) / (2.0f0 * step)
    @test ad.directional ≈ finite_difference rtol = 5.0f-2 atol = 1.0f-3
end
