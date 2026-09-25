using NeuralCrop
using Test

@testset "Atmospheric nitrogen deposition" begin
    soil = init_soil(2, soilparams.soildepth, identity)
    initial_nitrate = copy(soil.nitrogen.nitrate)
    initial_ammonium = copy(soil.nitrogen.ammonium)

    nitrogen_deposition!(
        soil,
        Float32[0.012, 0.034],
        Float32[0.056, 0.078],
    )

    @test soil.nitrogen.nitrate[1, :] == initial_nitrate[1, :] .+ Float32[0.012, 0.034]
    @test soil.nitrogen.ammonium[1, :] == initial_ammonium[1, :] .+ Float32[0.056, 0.078]
    @test soil.nitrogen.nitrate[2:end, :] == initial_nitrate[2:end, :]
    @test soil.nitrogen.ammonium[2:end, :] == initial_ammonium[2:end, :]
end
