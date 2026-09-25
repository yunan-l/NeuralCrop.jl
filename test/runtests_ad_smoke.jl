using NeuralCrop
using Enzyme
using Test

@testset "NeuralCrop Enzyme smoke" begin
    include("ad/test_adapter_contract.jl")
    include("ad/test_enzyme_adapter.jl")
end
