import CUDA
include(joinpath(@__DIR__, "run_global_wheat_cpu.jl"))

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: run_global_wheat_gpu.jl CONFIG_TOML")
    println(run_global_wheat(abspath(ARGS[1]); backend_override = :cuda))
end
