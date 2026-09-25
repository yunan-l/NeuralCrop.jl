@testset "Enzyme reverse seasonal gradient" begin
    template = _daily_transition_fixture()
    days = (11, 12, 13)
    values = _seasonal_production_values(template, days)
    context = ADSeasonContext(
        trues(3),
        (gpp = trues(3), reco = trues(3), et = trues(3)),
        (
            gpp = values.gpp .+ 0.1f0,
            reco = values.reco .+ 0.05f0,
            et = values.et .+ 0.05f0,
        ),
        (gpp = 1.0f0, reco = 1.0f0, et = 1.0f0),
    )
    parameter_names = (:gmin, :b)
    theta = Float32[cft1.gmin, cft1.b]
    direction = Float32[1, -1] ./ sqrt(2.0f0)
    primal = _seasonal_loss_value(template, theta, parameter_names, days, context)

    state = deepcopy(template.state)
    enzyme_prepare_daily_state!(state)
    gradient = zeros(Float32, length(theta))
    reverse_result = Enzyme.autodiff(
        Enzyme.ReverseWithPrimal,
        enzyme_seasonal_loss,
        Enzyme.Duplicated(theta, gradient),
        Enzyme.Duplicated(state, enzyme_zero_tangent(state)),
        Enzyme.Const(cft1),
        Enzyme.Const(template.global_parameters),
        Enzyme.Const(template.climate),
        Enzyme.Const(parameter_names),
        Enzyme.Const(days),
        Enzyme.Const(template.layer_depth),
        Enzyme.Const(context),
    )
    @test reverse_result[2] ≈ primal rtol = 1.0f-6 atol = 1.0f-7
    @test all(isfinite, gradient)

    forward = enzyme_forward_gradient(
        enzyme_seasonal_loss,
        theta,
        () -> begin
            fresh_state = deepcopy(template.state)
            return fresh_state, enzyme_zero_tangent(fresh_state)
        end,
        cft1,
        template.global_parameters,
        template.climate,
        parameter_names,
        days,
        template.layer_depth,
        context,
    )
    @test gradient ≈ forward rtol = 5.0f-2 atol = 1.0f-4
    @test dot(gradient, direction) ≈ dot(forward, direction) rtol = 1.0f-5 atol = 1.0f-6
end
