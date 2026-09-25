const _NeuralCropEnzymeExtension = Base.get_extension(NeuralCrop, :NeuralCropEnzymeExt)

function _model_parameter_replacement_sum(theta, base_parameters, parameter_names)
    parameters = _NeuralCropEnzymeExtension._replace_model_parameters(
        base_parameters, theta, parameter_names,
    )
    lpjml = parameters.lpjml
    return lpjml.k_soil10.fast + lpjml.PRIESTLEY_TAYLOR +
        lpjml.soildepth_evap + lpjml.soil_infil
end

@testset "Enzyme model-parameter replacement remains differentiable" begin
    parameter_names = (:k_soil10_fast, :PRIESTLEY_TAYLOR, :soildepth_evap, :soil_infil)
    base_parameters = ModelParameters(Float32)
    lpjml = base_parameters.lpjml
    theta = Float32[
        lpjml.k_soil10.fast,
        lpjml.PRIESTLEY_TAYLOR,
        lpjml.soildepth_evap,
        lpjml.soil_infil,
    ]
    gradient = zeros(Float32, length(theta))
    Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal),
        _model_parameter_replacement_sum,
        Enzyme.Duplicated(theta, gradient),
        Enzyme.Const(base_parameters),
        Enzyme.Const(parameter_names),
    )

    @test gradient == ones(Float32, length(theta))
end

@testset "Enzyme blockwise soil reverse" begin
    template = _daily_transition_fixture(14, 10)
    days = (11, 12, 13)
    values = _seasonal_production_values(template, days)
    context = ADSeasonContext(
        trues(length(days)),
        (gpp = falses(length(days)), reco = trues(length(days)), et = trues(length(days))),
        (gpp = values.gpp, reco = values.reco .+ 0.05f0, et = values.et .+ 0.05f0),
        (gpp = 1.0f0, reco = 1.0f0, et = 1.0f0),
    )
    parameter_names = (:k_soil10_fast, :PRIESTLEY_TAYLOR, :soildepth_evap, :soil_infil)
    parameters = ModelParameters(Float32)
    lpjml = parameters.lpjml
    theta = Float32[
        lpjml.k_soil10.fast,
        lpjml.PRIESTLEY_TAYLOR,
        lpjml.soildepth_evap,
        lpjml.soil_infil,
    ]
    state_factory() = begin
        state = deepcopy(template.state)
        return state, enzyme_zero_tangent(state)
    end

    result = enzyme_seasonal_soil_gradient_blockwise(
        theta,
        state_factory,
        cft1,
        parameters,
        parameter_names,
        template.climate,
        days,
        template.layer_depth,
        context;
        block_days = 2,
    )

    @test isfinite(result.primal)
    @test all(isfinite, result.gradient)
end

@testset "Enzyme joint CFT and soil reverse" begin
    template = _daily_transition_fixture(14, 10)
    days = (11, 12, 13)
    values = _seasonal_production_values(template, days)
    context = ADSeasonContext(
        trues(length(days)),
        (gpp = trues(length(days)), reco = trues(length(days)), et = trues(length(days))),
        (
            gpp = values.gpp .+ 0.05f0,
            reco = values.reco .+ 0.05f0,
            et = values.et .+ 0.05f0,
        ),
        (gpp = 1.0f0, reco = 1.0f0, et = 1.0f0),
    )
    cft_parameter_names = (:gmin, :b)
    soil_parameter_names = (
        :k_soil10_fast, :PRIESTLEY_TAYLOR, :soildepth_evap, :soil_infil,
    )
    theta_cft = Float32[cft1.gmin, cft1.b]
    parameters = ModelParameters(Float32)
    lpjml = parameters.lpjml
    theta_soil = Float32[
        lpjml.k_soil10.fast,
        lpjml.PRIESTLEY_TAYLOR,
        lpjml.soildepth_evap,
        lpjml.soil_infil,
    ]
    state_factory() = begin
        state = deepcopy(template.state)
        return state, enzyme_zero_tangent(state)
    end

    joint = enzyme_seasonal_joint_gradient_blockwise(
        theta_cft,
        theta_soil,
        state_factory,
        cft1,
        parameters,
        cft_parameter_names,
        soil_parameter_names,
        template.climate,
        days,
        template.layer_depth,
        context;
        block_days = 2,
    )
    cft_only = enzyme_seasonal_gradient_blockwise(
        theta_cft,
        state_factory,
        cft1,
        parameters,
        template.climate,
        cft_parameter_names,
        days,
        template.layer_depth,
        context;
        block_days = 2,
    )
    soil_only = enzyme_seasonal_soil_gradient_blockwise(
        theta_soil,
        state_factory,
        cft1,
        parameters,
        soil_parameter_names,
        template.climate,
        days,
        template.layer_depth,
        context;
        block_days = 2,
    )

    @test isfinite(joint.primal)
    @test all(isfinite, joint.cft_gradient)
    @test all(isfinite, joint.soil_gradient)
    @test joint.forward_primal ≈ cft_only.forward_primal rtol = 1.0f-6 atol = 1.0f-7
    @test joint.forward_primal ≈ soil_only.forward_primal rtol = 1.0f-6 atol = 1.0f-7
    @test joint.cft_gradient ≈ cft_only.gradient rtol = 1.0f-5 atol = 1.0f-6
    @test joint.soil_gradient ≈ soil_only.gradient rtol = 1.0f-5 atol = 1.0f-6
end
