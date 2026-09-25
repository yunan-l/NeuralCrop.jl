using NeuralCrop
using Test

@testset "AD season contract" begin
    days = 8
    context = ADSeasonContext(
        trues(days),
        (gpp = trues(days), reco = falses(days), et = trues(days)),
        (gpp = fill(2.0f0, days), reco = fill(NaN32, days), et = fill(1.0f0, days)),
        (gpp = 2.0f0, reco = 1.0f0, et = 1.0f0),
    )
    @test context.calendar_days == 365
    @test context.counts == (gpp = days, reco = 0, et = days)
    @test context.target_count == 2
    @test context.growth_mask == trues(days)

    @test_throws ArgumentError ADSeasonContext(
        trues(days),
        (gpp = trues(days), reco = trues(days), et = trues(days)),
        (gpp = fill(1.0f0, days), reco = fill(1.0f0, days), et = fill(1.0f0, days)),
        (gpp = 0.0f0, reco = 1.0f0, et = 1.0f0),
    )
    @test_throws ArgumentError ADSeasonContext(
        trues(days),
        (gpp = trues(days), reco = trues(days), et = trues(days)),
        (gpp = fill(1.0f0, days), reco = fill(1.0f0, days), et = fill(1.0f0, days)),
        (gpp = 1.0f0, reco = 1.0f0, et = 1.0f0);
        calendar_days = 366,
    )
end
