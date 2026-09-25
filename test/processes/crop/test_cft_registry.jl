@testset "LPJmL crop CFT registry" begin
    @test length(CFTS) == 12
    @test getfield.(CFTS, :name) == Tuple(1:12)
    @test CFT_NAMES[4] == "tropical cereals"
    @test CFT_NAMES[9] == "oil crops soybean"
    @test crop_cft(4) === cft4
    @test crop_cft("oil crops soybean") === cft9

    @test cft1.path == 1
    @test cft1.laimax == 7.0f0
    @test cft1.hlimit == 360
    @test cft4.path == 2
    @test cft4.fphusen == 0.85f0
    @test cft9.path == 1
    @test cft9.basetemp.low == 7.0f0
    @test cft12.path == 2
    @test cft12.hiopt == 0.8f0

    @test_throws ArgumentError crop_cft(13)
    @test_throws ArgumentError crop_cft("soybean")
end
