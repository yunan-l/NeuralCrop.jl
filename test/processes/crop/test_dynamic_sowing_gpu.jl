using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA dynamic sowing matches CPU for target CFT triggers" begin
    function run_trigger(device, cft, day; winter = false, temperature = (13f0, 15f0), wet = false)
        crop = init_crop(Float32, 1, device)
        climate = init_climbuf(Float32, 1, device)
        crop.auxiliary.calendar.prescribed_sowing_date .= Int32(100)
        crop.auxiliary.phenology.winter_type .= winter
        climate.mtemp20 .= device(reshape(Float32[2, 4, 15, 18, 20, 22, 23, 22, 18, 12, 7, 3], 12, 1))
        climate.mprec20 .= device(reshape(Float32[10, 10, 10, 10, 100, 100, 100, 100, 10, 10, 10, 10], 12, 1))
        climate.mpet20 .= device(fill(20f0, 12, 1))
        climate.temp[end - 1, :] .= temperature[1]
        climate.temp[end, :] .= temperature[2]
        wet && (climate.prec[day, :] .= 1f0)
        update_dynamic_sowing_calendar!(climate, cft, crop.auxiliary.phenology.winter_type)
        dynamic_sowing_date!(crop, climate, cft, day)
        return Array(crop.auxiliary.calendar.sowing_date)
    end

    for (cft, day, winter, temperature, wet) in (
        (cft1, 60, false, (4f0, 6f0), false),
        (cft3, 60, false, (13f0, 15f0), false),
        (cft2, 121, false, (13f0, 15f0), true),
        (cft9, 121, false, (13f0, 15f0), true),
    )
        @test run_trigger(CuArray, cft, day; winter, temperature, wet) ==
              run_trigger(identity, cft, day; winter, temperature, wet)
    end
end
