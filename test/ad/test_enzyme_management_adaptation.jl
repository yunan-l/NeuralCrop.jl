using LinearAlgebra

function _management_adaptation_state_factory()
    template = _daily_transition_fixture(18, 10)
    state = deepcopy(template.state)
    enzyme_prepare_daily_state!(state)
    return (; template..., state)
end

function _management_adaptation_loss(template, theta, context)
    state = deepcopy(template.state)
    enzyme_prepare_daily_state!(state)
    return enzyme_management_yield_loss(
        theta,
        state,
        cft1,
        template.global_parameters,
        template.climate,
        11:18,
        template.layer_depth,
        context,
    )
end

function _management_split_adaptation_loss(template, theta, context)
    state = deepcopy(template.state)
    enzyme_prepare_daily_state!(state)
    return enzyme_management_yield_split_loss(
        theta,
        state,
        cft1,
        template.global_parameters,
        template.climate,
        11:18,
        template.layer_depth,
        context,
    )
end

function _management_joint_adaptation_loss(template, theta, context)
    state = deepcopy(template.state)
    enzyme_prepare_daily_state!(state)
    return enzyme_joint_adaptation_yield_loss(
        theta,
        state,
        cft1,
        template.global_parameters,
        template.climate,
        11:18,
        template.layer_depth,
        context,
    )
end

@testset "Enzyme fixed-event management adaptation" begin
    template = _management_adaptation_state_factory()
    context = ManagementAdaptationContext((11, 14), (0.2f0, 0.8f0); nitrogen_cost = 0.05f0)
    theta = Float32[24.55]
    direction = Float32[1]
    primal = _management_adaptation_loss(template, theta, context)
    @test isfinite(primal)

    state_factory() = begin
        state = deepcopy(template.state)
        return state, enzyme_zero_tangent(state)
    end
    ad = enzyme_forward_directional(
        enzyme_management_yield_loss,
        theta,
        direction,
        state_factory,
        cft1,
        template.global_parameters,
        template.climate,
        11:18,
        template.layer_depth,
        context;
        return_primal = true,
    )
    @test ad.primal ≈ primal rtol = 1.0f-6 atol = 1.0f-7

    step = 1.0f-3
    finite_difference = (
        _management_adaptation_loss(template, theta .+ step .* direction, context) -
        _management_adaptation_loss(template, theta .- step .* direction, context)
    ) / (2.0f0 * step)
    @test ad.directional ≈ finite_difference rtol = 5.0f-2 atol = 1.0f-4

    @test_throws ArgumentError ManagementAdaptationContext((11, 11), (0.2f0, 0.8f0))
    @test_throws ArgumentError ManagementAdaptationContext((11, 14), (0.1f0, 0.8f0))
end


@testset "Enzyme joint management and cultivar adaptation" begin
    template = _management_adaptation_state_factory()
    context = ManagementAdaptationContext((11, 14), (0.2f0, 0.8f0); nitrogen_cost = 0.05f0)
    theta = Float32[24.55, 0.2, 1.0, 0.0, 1.0, 1.0]

    state_factory() = begin
        state = deepcopy(template.state)
        return state, enzyme_zero_tangent(state)
    end
    gradient = enzyme_forward_gradient(
        enzyme_joint_adaptation_yield_loss,
        theta,
        state_factory,
        cft1,
        template.global_parameters,
        template.climate,
        11:18,
        template.layer_depth,
        context,
    )
    @test all(isfinite, gradient)

    steps = Float32[1.0f-3, 1.0f-4, 1.0f-4, 1.0f-3, 1.0f-4, 1.0f-4]
    finite_difference = similar(theta)
    for index in eachindex(theta)
        perturbation = zeros(Float32, length(theta))
        perturbation[index] = steps[index]
        finite_difference[index] = (
            _management_joint_adaptation_loss(template, theta .+ perturbation, context) -
            _management_joint_adaptation_loss(template, theta .- perturbation, context)
        ) / (2.0f0 * steps[index])
    end
    @test all(isfinite, finite_difference)
    @test gradient ≈ finite_difference rtol = 5.0f-2 atol = 2.0f-4
end


@testset "Enzyme two-event fertilizer allocation" begin
    template = _management_adaptation_state_factory()
    context = ManagementAdaptationContext((11, 14), (0.2f0, 0.8f0); nitrogen_cost = 0.05f0)
    theta = Float32[24.55, 0.2]

    state_factory() = begin
        state = deepcopy(template.state)
        return state, enzyme_zero_tangent(state)
    end
    gradient = enzyme_forward_gradient(
        enzyme_management_yield_split_loss,
        theta,
        state_factory,
        cft1,
        template.global_parameters,
        template.climate,
        11:18,
        template.layer_depth,
        context,
    )
    @test all(isfinite, gradient)

    for index in eachindex(theta)
        step = index == 1 ? 1.0f-3 : 1.0f-4
        direction = zeros(Float32, length(theta))
        direction[index] = step
        finite_difference = (
            _management_split_adaptation_loss(template, theta .+ direction, context) -
            _management_split_adaptation_loss(template, theta .- direction, context)
        ) / (2.0f0 * step)
        @test gradient[index] ≈ finite_difference rtol = 5.0f-2 atol = 2.0f-4
    end

    three_event_context = ManagementAdaptationContext(
        (11, 13, 15), (0.2f0, 0.3f0, 0.5f0),
    )
    @test_throws ArgumentError _management_split_adaptation_loss(
        template, theta, three_event_context,
    )
end
