function _irrigated_daily_transition_objective(
    theta,
    state,
    base_cft,
    global_parameters,
    climate,
    parameter_names,
    day,
    observable,
    layer_depth,
)
    return enzyme_daily_transition_objective(
        theta,
        state,
        base_cft,
        global_parameters,
        climate,
        parameter_names,
        day,
        observable,
        layer_depth;
        irrigation = true,
    )
end

function _set_frozen_irrigated_state!(state)
    water = NeuralCrop.soil_water_prognostic(state)
    field_capacity = NeuralCrop.soil_water_auxiliary(state).field_capacity
    layer_depth = NeuralCrop.soil_properties(state).layer_depth
    target_storage = field_capacity .* reshape(layer_depth, :, 1)
    water.storage .= target_storage .* 0.5f0
    water.ice_storage .= target_storage .* 0.25f0
    NeuralCrop.partition_soil_water_ice!(state)
    return target_storage
end

@testset "Enzyme irrigated frozen-soil reset" begin
    template = _daily_transition_fixture()
    extension = Base.get_extension(NeuralCrop, :NeuralCropEnzymeExt)
    @test extension !== nothing

    production_state = deepcopy(template.state)
    enzyme_state = deepcopy(template.state)
    target_storage = _set_frozen_irrigated_state!(production_state)
    _set_frozen_irrigated_state!(enzyme_state)

    soil_evapotranspiration!(production_state, production_state; irrigation = true)
    extension._enzyme_soil_evapotranspiration!(enzyme_state; irrigation = true)

    production_water = NeuralCrop.soil_water_prognostic(production_state)
    enzyme_water = NeuralCrop.soil_water_prognostic(enzyme_state)
    @test enzyme_water.storage == production_water.storage
    @test enzyme_water.ice_storage == production_water.ice_storage
    @test enzyme_water.storage + enzyme_water.ice_storage ≈ target_storage

    parameter_names = (:gmin, :b)
    theta = Float32[cft1.gmin, cft1.b]
    state = deepcopy(template.state)
    enzyme_prepare_daily_state!(state)
    target_storage = _set_frozen_irrigated_state!(state)
    primal = _irrigated_daily_transition_objective(
        theta,
        state,
        cft1,
        template.global_parameters,
        template.climate,
        parameter_names,
        template.day,
        :et,
        template.layer_depth,
    )
    water = NeuralCrop.soil_water_prognostic(state)
    @test isfinite(primal)
    @test water.storage + water.ice_storage ≈ target_storage

    reverse_state = deepcopy(template.state)
    enzyme_prepare_daily_state!(reverse_state)
    _set_frozen_irrigated_state!(reverse_state)
    gradient = zeros(Float32, length(theta))
    reverse_result = Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal),
        _irrigated_daily_transition_objective,
        Enzyme.Duplicated(theta, gradient),
        Enzyme.Duplicated(reverse_state, enzyme_zero_tangent(reverse_state)),
        Enzyme.Const(cft1),
        Enzyme.Const(template.global_parameters),
        Enzyme.Const(template.climate),
        Enzyme.Const(parameter_names),
        Enzyme.Const(template.day),
        Enzyme.Const(:et),
        Enzyme.Const(template.layer_depth),
    )
    @test isfinite(reverse_result[2])
    @test all(isfinite, gradient)

    forward_gradient = enzyme_forward_gradient(
        _irrigated_daily_transition_objective,
        theta,
        () -> begin
            fresh_state = deepcopy(template.state)
            enzyme_prepare_daily_state!(fresh_state)
            _set_frozen_irrigated_state!(fresh_state)
            return fresh_state, enzyme_zero_tangent(fresh_state)
        end,
        cft1,
        template.global_parameters,
        template.climate,
        parameter_names,
        template.day,
        :et,
        template.layer_depth,
    )
    @test gradient ≈ forward_gradient rtol = 5.0f-2 atol = 1.0f-4
end
