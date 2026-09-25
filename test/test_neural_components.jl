@testset "NeuralCrop component contract" begin
    components = NeuralComponents()

    @test components.gpp
    @test !components.gpp_residual
    @test !components.lambda
    @test !components.vcmax
    @test components.respiration
    @test components.transpiration
    @test !components.allocation
    @test !components.decomposition
    @test !components.snowmelt
    @test !components.evaporation

    layout = NeuralCropLayout(; components)
    @test neural_parameter_count(layout) == 42_127
    @test neural_trainable_parameter_count(layout) == 14_275

    theta = initialize_neural_parameters(layout; T = Float32, seed = 20260921)
    gpp = neural_gpp(
        theta, layout, 12.0f0, 1.0f7, 0.7f0, 3.0f0, 4.0f0,
        20.0f0, 410.0f0, 0.6f0,
    )
    @test 0.0f0 < gpp < 50.0f0
    @test neural_crop_respiration_multiplier(
        theta, layout, 10.0f0, 1.0f0, 20.0f0, 15.0f0,
        100.0f0, 50.0f0, 20.0f0,
    ) == 1.0f0
    @test neural_transpiration_multiplier(
        theta, layout, 12.0f0, 4.0f0, 0.7f0, 0.2f0,
        0.5f0, 100.0f0, 0.6f0,
    ) == 1.0f0

    residual_components = NeuralComponents(; gpp_residual = true)
    residual_layout = NeuralCropLayout(; components = residual_components)
    residual_theta = initialize_neural_parameters(
        residual_layout; T = Float32, seed = 20260921,
    )
    multiplier = neural_gpp_multiplier(
        residual_theta, residual_layout, 12.0f0, 1.0f7, 0.7f0, 3.0f0,
        4.0f0, 20.0f0, 410.0f0, 0.6f0,
    )
    @test multiplier == 1.0f0
    @test 0.0f0 <= multiplier <= 2.0f0
    @test_throws ArgumentError NeuralCropLayout(;
        components = NeuralComponents(; gpp = false, gpp_residual = true),
    )
end
