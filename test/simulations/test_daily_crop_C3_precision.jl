using NeuralCrop
using Test

precision_backend(::typeof(identity), values) = values
precision_backend(device, values) = device(values)

function run_c3_precision_smoke(::Type{T}, device = identity) where {T <: AbstractFloat}
    cells = 1
    layers = 5
    days = 10
    backend(values) = precision_backend(device, values)
    initial_data = (
        latitude = backend(T[45]),
        soilparams = (
            ph = backend(T[6.5]),
            w_sat = backend(fill(T(0.45), layers, cells)),
            sand = backend(reshape(T[0.4], 1, cells)),
            clay = backend(reshape(T[0.2], 1, cells)),
            tdiff_0 = backend(T[0.7]),
            tdiff_15 = backend(T[0.75]),
            soildepth = T[200, 300, 500, 1000, 1000],
        ),
        ModelState = (
            crop = (
                sdate = backend(Int32[1]),
                phu = backend(T[543]),
                manure = backend(zeros(T, cells)),
                fertilizer = backend(T[24.55]),
                residuefrac = backend(T[0.67]),
            ),
            u0 = (
                swc = backend(reshape(T[57.41, 55.32, 126.13, 274.59, 285.71], layers, cells)),
                litc = backend(reshape(T[0.13, 187.5, 225.36], 3, cells)),
                fastc = backend(reshape(T[548.97, 368.27, 313.79, 377.55, 344.65], layers, cells)),
                slowc = backend(reshape(T[1218.62, 753.33, 660.10, 792.63, 736.38], layers, cells)),
                litn = backend(reshape(T[0.0047, 6.47, 9.47], 3, cells)),
                fastn = backend(reshape(T[36.60, 24.55, 20.92, 25.17, 22.98], layers, cells)),
                slown = backend(reshape(T[81.24, 50.22, 44.01, 52.84, 49.09], layers, cells)),
            ),
        ),
    )

    climbuf, crop, pet, soil, managed_land, weather, output = init_states!(
        cft1, initial_data, cells, device; T = T,
    )
    climbuf.atemp .= T(10)
    climbuf.temp .= T(10)
    climbuf.atemp_mean .= T(10)
    climate = (
        temp = backend(fill(T(15), days, cells)),
        prec = backend(fill(T(1), days, cells)),
        sw = backend(fill(T(180), days, cells)),
        lw = backend(fill(T(-40), days, cells)),
        wind = backend(fill(T(2), days, cells)),
        co2 = backend(T[400]),
    )

    state = model_state(climbuf, crop, pet, soil, managed_land, weather, output)
    processes = ProcessModules(convert_precision(T, cft1), ModelParameters(T))
    daily_crop_C3!(
        1, days, processes, climate, state;
        fertilizer = :yes,
        nitrogen_limit_vcmax = false,
    )
    return (; crop, soil, output)
end

@testset "C3 selectable Float32/Float64 CPU precision" begin
    result32 = run_c3_precision_smoke(Float32)
    result64 = run_c3_precision_smoke(Float64)

    @test eltype(result32.output.crop.npp) == Float32
    @test eltype(result64.output.crop.npp) == Float64
    @test all(isfinite, result32.output.crop.npp)
    @test all(isfinite, result64.output.crop.npp)
    @test all(0.0f0 .<= result32.output.crop.water_deficit .<= 100.0f0)
    @test all(0.0 .<= result64.output.crop.water_deficit .<= 100.0)
    @test minimum(result32.output.crop.water_deficit) < 100.0f0
    @test result32.output.crop.water_deficit[end, :] ==
        result32.crop.auxiliary.stress.water_deficit
    @test result64.output.crop.water_deficit[end, :] ==
        result64.crop.auxiliary.stress.water_deficit
    @test all(isfinite, result64.soil.thermal.temperature)
    @test eltype(result64.soil.management.tillage_density_factor) == Float64
    # This compares two numerical precisions, not CPU/GPU execution at one
    # precision. Lambda bisection and nonlinear photosynthesis can amplify
    # rounding differences over successive days, so require the same physical
    # trajectory without demanding bit-level agreement across Float32/Float64.
    @test Float64.(result32.output.crop.npp) ≈ result64.output.crop.npp rtol = 1e-2 atol = 1e-5
    @test Float64.(result32.soil.water.storage) ≈ result64.soil.water.storage rtol = 5e-4 atol = 1e-5
end
