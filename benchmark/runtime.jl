using NeuralCrop
using BenchmarkTools
using JLD2

example_dir = joinpath(@__DIR__, "..", "examples")
initial_data = load(joinpath(example_dir, "initial_wheat.jld2"), "initial_data")
raw_climate = load(joinpath(example_dir, "climate_2000_2009.jld2"), "climate")
cells = length(initial_data.latitude)
days = min(365, size(raw_climate.temp, 1))
indices = collect(1:cells)
climate = ClimateDataLoader(raw_climate, indices, identity; T = Float32)

function benchmark_simulation()
    simulation = initialize_simulation(
        cft1, initial_data;
        indices,
        T = Float32,
        days,
        diagnostics = false,
        fertilizer = :yes,
    )
    run_simulation!(simulation, climate; end_day = days, spinup = false)
    return simulation
end

# Compile before measuring steady-state runtime.
benchmark_simulation()
trial = @benchmark benchmark_simulation() samples = 5 evals = 1
median_seconds = BenchmarkTools.median(trial).time / 1.0e9
println("cells = ", cells)
println("days = ", days)
println("median_seconds = ", median_seconds)
println("cell_days_per_second = ", cells * days / median_seconds)
println("median_allocated_bytes = ", BenchmarkTools.median(trial).memory)
println("projected_global_retained_memory = ", estimate_memory(
    67_420, 365;
    T = Float32,
    diagnostics = false,
    block_days = 31,
    backend = :accelerator,
))
global_stream = OutputStream(
    [
        OutputVariable(:crop, :npp),
        OutputVariable(:crop, :lai),
        OutputVariable(:crop, :yield),
    ];
    frequency = :monthly,
    cell_ids = 1:67_420,
)
println("projected_global_streaming_memory = ", estimate_memory(
    67_420, 365;
    T = Float32,
    diagnostics = false,
    block_days = 31,
    backend = :accelerator,
    prefetch = true,
    output_stream = global_stream,
))
