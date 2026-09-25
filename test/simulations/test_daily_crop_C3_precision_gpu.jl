using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

# Defines the shared deterministic CPU/GPU precision fixture and also runs its
# CPU-only regression before the GPU comparisons below.
include("test_daily_crop_C3_precision.jl")

function test_pointwise_approx(actual, expected; label, rtol, atol)
    absolute_error = abs.(actual .- expected)
    tolerance = atol .+ rtol .* max.(abs.(actual), abs.(expected))
    tolerance_ratio = absolute_error ./ max.(tolerance, eps(eltype(actual)))
    worst_index = argmax(tolerance_ratio)

    @info "CPU-GPU pointwise comparison" label eltype = eltype(actual) worst_index actual = actual[worst_index] expected = expected[worst_index] maximum_absolute_error = absolute_error[worst_index] maximum_tolerance_ratio = tolerance_ratio[worst_index] rtol atol

    @test all(absolute_error .<= tolerance)
end

@testset "CUDA C3 Float32/Float64 precision support" begin
    for T in (Float32, Float64)
        cpu = run_c3_precision_smoke(T, identity)
        gpu = run_c3_precision_smoke(T, CuArray)

        gpu_npp = Array(gpu.output.crop.npp)
        gpu_water = Array(gpu.soil.water.storage)
        @test eltype(gpu.output.crop.npp) == T
        @test eltype(gpu.soil.water.storage) == T
        @test all(isfinite, gpu_npp)
        @test all(isfinite, gpu_water)

        # Float32 GPU arithmetic may use fused operations and device math whose
        # rounding differs slightly from the CPU. Photosynthesis then propagates
        # that difference through the daily nonlinear state update. Keep a
        # strict, explicit bound instead of requiring bitwise-equivalent paths.
        npp_rtol = T === Float32 ? T(5e-4) : T(5e-5)
        test_pointwise_approx(
            gpu_npp,
            cpu.output.crop.npp;
            label = "C3 NPP",
            rtol = npp_rtol,
            atol = T(5e-6),
        )
        test_pointwise_approx(
            gpu_water,
            cpu.soil.water.storage;
            label = "soil water storage",
            rtol = T(5e-5),
            atol = T(5e-6),
        )
    end
end
