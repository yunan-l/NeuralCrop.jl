using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "Backend-neutral global inputs initialize and run on CUDA" begin
    T = Float32
    layers = 5
    initial = (
        backend_neutral = true,
        coords = Int32[1],
        latitude = T[45],
        crop = (
            sdate = Int32[1],
            phu = T[543],
            manure = T[0],
            fertilizer = T[20],
            residuefrac = T[0.67],
        ),
        soilparam = (
            soilph = T[6.5],
            w_sat = fill(T(0.45), layers, 1),
            sand = T[0.4],
            clay = T[0.2],
            tdiff_0 = T[0.7],
            tdiff_15 = T[0.75],
            soildepth = T[200, 300, 500, 1000, 1000],
        ),
        initial_state = (
            swc = reshape(T[57, 55, 126, 275, 286], layers, 1),
            litc = reshape(T[0.1, 188, 225], 3, 1),
            fastc = fill(T(10), layers, 1),
            slowc = fill(T(100), layers, 1),
            litn = reshape(T[0.005, 6.5, 9.5], 3, 1),
            fastn = fill(T(1), layers, 1),
            slown = fill(T(10), layers, 1),
        ),
    )
    simulation = initialize_simulation(
        cft1, initial;
        device = CuArray, T, days = 1, diagnostics = false, fertilizer = :no,
    )
    @test simulation.state.prognostic.soil.carbon.fast isa CuArray{T, 2}
    @test simulation.state.inputs.crop.calendar.sowing_date isa CuArray{Int32, 1}
    @test simulation.state.inputs.crop.calendar.prescribed_sowing_date isa CuArray{Int32, 1}
    @test simulation.managed_land.latitude isa CuArray{T, 1}

    climate = (
        backend_neutral = true,
        temp = fill(T(15), 1, 1),
        prec = fill(T(1), 1, 1),
        sw = fill(T(180), 1, 1),
        lw = fill(T(-40), 1, 1),
        co2 = T[400],
    )
    management = (
        sdate = Int32[1],
        phu = T[-700],
        manure = T[0],
        fertilizer = T[0],
        residuefrac = T[0.5],
    )
    run_simulation!(simulation, climate; spinup = false, management)
    @test simulation.simulated_days == 1
    @test Array(simulation.state.inputs.crop.phenology.phu) == T[700]
    @test Array(simulation.state.inputs.crop.phenology.winter_type) == Bool[true]
end
