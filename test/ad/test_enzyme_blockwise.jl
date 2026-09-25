@testset "Enzyme blockwise 365-day reverse" begin
    template = _daily_transition_fixture(375, 10)
    days = Tuple(11:375)
    _, context = _fixed_365_context(template, days)
    parameter_names = (:gmin, :b)
    theta = Float32[cft1.gmin, cft1.b]
    direction = Float32[1, -1] ./ sqrt(2.0f0)

    full_state = deepcopy(template.state)
    enzyme_prepare_daily_state!(full_state)
    full_gradient = zeros(Float32, length(theta))
    full_result = Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal),
        enzyme_seasonal_loss,
        Enzyme.Duplicated(theta, full_gradient),
        Enzyme.Duplicated(full_state, enzyme_zero_tangent(full_state)),
        Enzyme.Const(cft1),
        Enzyme.Const(template.global_parameters),
        Enzyme.Const(template.climate),
        Enzyme.Const(parameter_names),
        Enzyme.Const(days),
        Enzyme.Const(template.layer_depth),
        Enzyme.Const(context),
    )

    state_factory() = begin
        state = deepcopy(template.state)
        enzyme_prepare_daily_state!(state)
        return state, enzyme_zero_tangent(state)
    end
    blockwise = enzyme_seasonal_gradient_blockwise(
        theta,
        state_factory,
        cft1,
        template.global_parameters,
        template.climate,
        parameter_names,
        days,
        template.layer_depth,
        context;
        block_days = 30,
    )

    relative_error = abs.(blockwise.gradient .- full_gradient) ./
        max.(abs.(full_gradient), 1.0f-6)
    max_relative_error = maximum(relative_error)
    @info "365-day full versus blockwise reverse" full_gradient = Tuple(full_gradient) blockwise_gradient = Tuple(blockwise.gradient) relative_error = Tuple(relative_error) max_relative_error = max_relative_error full_loss = full_result[2] blockwise_loss = blockwise.primal blocks = length(blockwise.block_ranges)

    @test length(blockwise.block_ranges) == 13
    @test length(blockwise.block_ranges[end]) == 5
    @test blockwise.forward_primal ≈ full_result[2] rtol = 1.0f-6 atol = 1.0f-7
    @test blockwise.primal ≈ full_result[2] rtol = 2.0f-5 atol = 2.0f-6
    @test blockwise.gradient ≈ full_gradient rtol = 2.0f-3 atol = 2.0f-5
    @test dot(blockwise.gradient, direction) ≈ dot(full_gradient, direction) rtol = 2.0f-3 atol = 2.0f-5
end
