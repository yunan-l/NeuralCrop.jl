using Test

gpu_tests = String[]
for (directory, _, files) in walkdir(@__DIR__)
    for file in files
        endswith(file, "_gpu.jl") || continue
        file == basename(@__FILE__) && continue
        push!(gpu_tests, joinpath(directory, file))
    end
end
sort!(gpu_tests)

@testset "NeuralCrop GPU" begin
    for path in gpu_tests
        @testset "$(relpath(path, @__DIR__))" begin
            # Isolate helpers and constants declared by standalone GPU tests.
            test_module = Module(gensym(:NeuralCropGPUTest), true, true)
            Core.eval(
                test_module,
                :(include(path::AbstractString) = Base.include(@__MODULE__, path)),
            )
            Core.eval(test_module, :(using NeuralCrop))
            Base.include(test_module, joinpath(@__DIR__, "helpers", "model_state_fixture.jl"))
            Base.include(test_module, path)
        end
    end
end
