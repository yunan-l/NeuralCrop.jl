using NeuralCrop
using Test

@testset "Optional daily wind forcing" begin
    weather = init_weather(2, identity)
    @test weather.wind == fill(lpjmlparams.volatil_wind, 2)

    climate_without_wind = (
        temp = Float32[10 11; 12 13],
        prec = zeros(Float32, 2, 2),
        sw = zeros(Float32, 2, 2),
        lw = zeros(Float32, 2, 2),
        co2 = Float32[400],
    )
    readclimate!(climate_without_wind, weather, 1)
    @test weather.wind == fill(lpjmlparams.volatil_wind, 2)
    @test weather.no3_deposition == zeros(Float32, 2)
    @test weather.nh4_deposition == zeros(Float32, 2)

    climate_with_wind = merge(
        climate_without_wind,
        (wind = Float32[2.0 3.0; 4.0 5.0],),
    )
    readclimate!(climate_with_wind, weather, 2)
    @test weather.wind == Float32[4, 5]

    # Moving back to a legacy archive must not retain the previous day's wind.
    readclimate!(climate_without_wind, weather, 2)
    @test weather.wind == fill(lpjmlparams.volatil_wind, 2)

    climate_with_deposition = merge(
        climate_without_wind,
        (
            no3_deposition = Float32[0.001 0.002; 0.003 0.004],
            nh4_deposition = Float32[0.005 0.006; 0.007 0.008],
        ),
    )
    readclimate!(climate_with_deposition, weather, 2)
    @test weather.no3_deposition == Float32[0.003, 0.004]
    @test weather.nh4_deposition == Float32[0.007, 0.008]

    loader_input = (
        temp_spinup = zeros(Float32, 2, 2),
        temp = zeros(Float32, 2, 2),
        prec = zeros(Float32, 2, 2),
        swdown = zeros(Float32, 2, 2),
        lwnet = zeros(Float32, 2, 2),
        co2 = Float32[400],
        windspeed = Float32[1 2; 3 4],
        no3_deposition = Float32[0.001 0.002; 0.003 0.004],
        nh4_deposition = Float32[0.005 0.006; 0.007 0.008],
    )
    loaded = ClimateDataLoader(loader_input, [2], identity)
    @test hasproperty(loaded, :wind)
    @test loaded.wind == Float32[2; 4;;]
    @test loaded.no3_deposition == Float32[0.002; 0.004;;]
    @test loaded.nh4_deposition == Float32[0.006; 0.008;;]

    loaded64 = ClimateDataLoader(loader_input, [2], identity; T = Float64)
    @test eltype(loaded64.temp) == Float64
    @test eltype(loaded64.co2) == Float64
    @test eltype(loaded64.wind) == Float64
    @test eltype(loaded64.no3_deposition) == Float64

    monthly_deposition = reshape(Float32.(1:12), 12, 1)
    monthly_loader_input = merge(loader_input, (
        no3_deposition = monthly_deposition,
        nh4_deposition = monthly_deposition .* 2,
        temp_spinup = zeros(Float32, 365, 2),
        temp = zeros(Float32, 365, 2),
        prec = zeros(Float32, 365, 2),
        swdown = zeros(Float32, 365, 2),
        lwnet = zeros(Float32, 365, 2),
    ))
    monthly_loaded = ClimateDataLoader(monthly_loader_input, [1], identity)
    @test size(monthly_loaded.no3_deposition) == (365, 1)
    # LPJmL's 1 January lies 17/31 of the way from December to January.
    @test monthly_loaded.no3_deposition[1, 1] ≈ 12.0f0 + 17.0f0 / 31.0f0 * (1.0f0 - 12.0f0)
    # The mid-month value is exactly the current monthly value.
    @test monthly_loaded.no3_deposition[15, 1] == 1.0f0
    @test monthly_loaded.nh4_deposition ≈ monthly_loaded.no3_deposition .* 2
end


@testset "Climate read kernel is deterministic across independent state" begin
    cells = 5
    days = 3
    climate = (
        temp = reshape(Float32.(1:(days * cells)), days, cells),
        prec = reshape(Float32.(11:(10 + days * cells)), days, cells),
        sw = reshape(Float32.(101:(100 + days * cells)), days, cells),
        lw = reshape(Float32.(-20:(-21 + days * cells)), days, cells),
        wind = reshape(Float32.(range(1, 6; length = days * cells)), days, cells),
        co2 = Float32[400, 405],
    )
    reference = init_weather(cells, identity)
    kernel = init_weather(cells, identity)
    reference_co2 = NeuralCrop.readclimate!(climate, reference, 2)
    kernel_co2 = readclimate!(climate, kernel, 2)
    for field in (:temp, :prec, :swr, :lwr, :wind, :no3_deposition, :nh4_deposition, :annual_co2)
        @test getproperty(kernel, field) ≈ getproperty(reference, field)
    end
    @test reference_co2 === reference.annual_co2
    @test kernel_co2 === kernel.annual_co2

    daily_co2 = reshape(Float32.(range(390, 430; length = days * cells)), days, cells)
    daily_climate = merge(climate, (co2 = daily_co2,))
    reference_co2 = NeuralCrop.readclimate!(daily_climate, reference, 3)
    kernel_co2 = readclimate!(daily_climate, kernel, 3)
    for field in (:temp, :prec, :swr, :lwr, :wind, :no3_deposition, :nh4_deposition, :daily_co2)
        @test getproperty(kernel, field) ≈ getproperty(reference, field)
    end
    @test reference_co2 === reference.daily_co2
    @test kernel_co2 === kernel.daily_co2

    global_daily_climate = merge(climate, (
        co2 = Float32[401, 402, 403],
        co2_daily = true,
    ))
    reference_co2 = NeuralCrop.readclimate!(global_daily_climate, reference, 2)
    kernel_co2 = readclimate!(global_daily_climate, kernel, 2)
    @test reference_co2[1] == 40.2f0
    @test kernel_co2 == reference_co2
end
