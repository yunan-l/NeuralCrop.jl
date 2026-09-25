using NeuralCrop
using Test

@testset "Annual climate diagnostics match column-wise reference" begin
    month_lengths = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    monthly = Float32[
        8  -4  12;
        6  -2  10;
        4   0   8;
        2   2   6;
        0   4   4;
       -2   6   2;
       -4   8   0;
       -6  10  -2;
       -8  12  -4;
      -10  14  -6;
      -12  16  -8;
      -14  18 -10;
    ]
    daily = zeros(Float32, 365, 3)
    first_day = 1
    for month in 1:12
        last_day = first_day + month_lengths[month] - 1
        @views daily[first_day:last_day, :] .= reshape(monthly[month, :], 1, :)
        first_day = last_day + 1
    end

    climbuf = init_climbuf(3, identity)
    annual_climbuf!(daily, climbuf, cft1)

    expected_minimum = sort(monthly; dims = 1)[1:5, :]
    expected_annual_mean = vec(sum(
        monthly .* reshape(collect(Float32, month_lengths), :, 1); dims = 1,
    )) ./ 365.0f0
    expected_vreq = zeros(Float32, 3)
    for cell in 1:3, rank in 1:5
        temperature = expected_minimum[rank, cell]
        if temperature <= cft1.tv_opt.low && temperature > -9999
            expected_vreq[cell] += cft1.pvd_max / 5
        elseif temperature < cft1.tv_opt.high
            expected_vreq[cell] += cft1.pvd_max / 5 *
                (1 - (temperature - cft1.tv_opt.low) /
                     (cft1.tv_opt.high - cft1.tv_opt.low))
        end
    end

    @test climbuf.mtemp ≈ monthly
    @test climbuf.mtemp20 ≈ monthly
    @test climbuf.min_temp ≈ expected_minimum
    @test climbuf.V_req_a ≈ expected_vreq
    @test climbuf.V_req ≈ expected_vreq
    @test climbuf.atemp_mean ≈ expected_annual_mean
end

@testset "Fixed prescribed management freezes V_req but not climate history" begin
    daily = fill(5.0f0, 365, 1)
    climbuf = init_climbuf(1, identity)
    annual_climbuf!(daily, climbuf, cft1)
    vreq = copy(climbuf.V_req)
    vreq_a = copy(climbuf.V_req_a)

    annual_climbuf!(
        fill(15.0f0, 365, 1), climbuf, cft1;
        update_vernalization_requirement = false,
    )

    @test climbuf.V_req == vreq
    @test climbuf.V_req_a == vreq_a
    @test climbuf.mtemp20 != fill(5.0f0, 12, 1)
end
