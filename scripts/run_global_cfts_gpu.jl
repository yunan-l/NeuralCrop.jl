import CUDA
include(joinpath(@__DIR__, "run_global_cfts_cpu.jl"))

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (1, 3) || error(
        "usage: run_global_cfts_gpu.jl CONFIG_TOML [CFT_ID rainfed|irrigated]",
    )
    cft_id = length(ARGS) == 3 ? parse(Int, ARGS[2]) : nothing
    water = length(ARGS) == 3 ? Symbol(lowercase(ARGS[3])) : nothing
    length(ARGS) == 3 && water in (:rainfed, :irrigated) ||
        length(ARGS) == 1 || error("water system must be rainfed or irrigated")
    println(run_global_cfts(
        abspath(ARGS[1]);
        backend_override = :cuda,
        cft_id,
        irrigated = isnothing(water) ? nothing : water === :irrigated,
    ))
end
