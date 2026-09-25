using NCDatasets

function _write_daily_climate(path, variable_name, units, values; year = 2015)
    NCDataset(path, "c") do ds
        defDim(ds, "time", 365)
        defDim(ds, "latitude", 1)
        defDim(ds, "longitude", 1)
        time = defVar(ds, "time", Int32, ("time",))
        time[:] = Int32.(0:364)
        time.attrib["units"] = "days since $year-01-01"
        time.attrib["calendar"] = "noleap"
        latitude = defVar(ds, "latitude", Float64, ("latitude",))
        longitude = defVar(ds, "longitude", Float64, ("longitude",))
        latitude[:] = [0.5]
        longitude[:] = [0.5]
        variable = defVar(ds, variable_name, Float32, ("time", "latitude", "longitude"))
        variable[:, 1, 1] = values
        variable.attrib["units"] = units
    end
    return path
end

function _deposition_catalog(directory, climate_year, deposition_values)
    temp = _write_daily_climate(
        joinpath(directory, "temp.nc"), "temp", "degC", fill(10.0f0, 365);
        year = climate_year,
    )
    prec = _write_daily_climate(
        joinpath(directory, "prec.nc"), "prec", "mm/day", zeros(Float32, 365);
        year = climate_year,
    )
    lwnet = _write_daily_climate(
        joinpath(directory, "lwnet.nc"), "lwnet", "W/m2", fill(-50.0f0, 365);
        year = climate_year,
    )
    swdown = _write_daily_climate(
        joinpath(directory, "swdown.nc"), "swdown", "W/m2", fill(100.0f0, 365);
        year = climate_year,
    )
    no3 = _write_monthly_deposition(
        joinpath(directory, "no3.nc"), "noy", deposition_values,
    )
    nh4 = _write_monthly_deposition(
        joinpath(directory, "nh4.nc"), "nhx", 2.0f0 .* deposition_values,
    )
    co2 = joinpath(directory, "co2.txt")
    write(co2, "$climate_year 400.0\n")
    catalog = DatasetCatalog(Dict(
        :temp => DatasetSpec(temp, "temp"; units = "degC"),
        :prec => DatasetSpec(prec, "prec"; units = "mm/day"),
        :lwnet => DatasetSpec(lwnet, "lwnet"; units = "W/m2"),
        :swdown => DatasetSpec(swdown, "swdown"; units = "W/m2"),
        :co2 => DatasetSpec(co2, "co2"; units = "ppm"),
        :no3_deposition => DatasetSpec(no3, "noy"; units = "g/m2/day"),
        :nh4_deposition => DatasetSpec(nh4, "nhx"; units = "g/m2/day"),
    ), CFTRegistry([1], ["test_crop"]))
    return catalog, co2
end

function _write_monthly_deposition(path, variable_name, values)
    NCDataset(path, "c") do ds
        defDim(ds, "time", 12)
        defDim(ds, "latitude", 1)
        defDim(ds, "longitude", 1)
        time = defVar(ds, "time", Int32, ("time",))
        time[:] = Int32[30, 58, 89, 119, 150, 180, 211, 242, 272, 303, 333, 364]
        time.attrib["units"] = "days since 2015-01-01"
        time.attrib["calendar"] = "noleap"
        latitude = defVar(ds, "latitude", Float64, ("latitude",))
        longitude = defVar(ds, "longitude", Float64, ("longitude",))
        latitude[:] = [0.5]
        longitude[:] = [0.5]
        variable = defVar(ds, variable_name, Float32, ("time", "latitude", "longitude"))
        variable[:, 1, 1] = values
        variable.attrib["units"] = "g/m2/day"
    end
    return path
end

@testset "Monthly nitrogen deposition climate forcing" begin
    mktempdir() do directory
        catalog, co2 = _deposition_catalog(directory, 2015, Float32.(1:12))
        grid = GridIndex([0.5], [0.5], reshape(Int32[0], 1, 1), Int32[0], Int32[1], Int32[1])
        block = only(collect(climate_blocks(catalog, grid; co2_path = co2, block_days = 365)))
        forcing = climate_forcing(block)

        @test size(forcing.no3_deposition) == (365, 1)
        @test size(forcing.nh4_deposition) == (365, 1)
        @test forcing.no3_deposition[15, 1] == 1.0f0
        @test forcing.nh4_deposition[15, 1] == 2.0f0
        @test forcing.no3_deposition[1, 1] ≈ 12.0f0 + (17.0f0 / 31.0f0) * (1.0f0 - 12.0f0)
        @test forcing.nh4_deposition[1, 1] ≈ 2.0f0 * forcing.no3_deposition[1, 1]
    end
end

@testset "Fixed-year nitrogen deposition forcing" begin
    mktempdir() do directory
        catalog, co2 = _deposition_catalog(directory, 2016, Float32.(1:12))
        grid = GridIndex([0.5], [0.5], reshape(Int32[0], 1, 1), Int32[0], Int32[1], Int32[1])
        block = only(collect(climate_blocks(
            catalog, grid;
            co2_path = co2,
            block_days = 365,
            nitrogen_deposition_year = 2015,
        )))
        forcing = climate_forcing(block)

        @test forcing.no3_deposition[15, 1] == 1.0f0
        @test forcing.no3_deposition[365, 1] ≈
            12.0f0 + (16.0f0 / 31.0f0) * (1.0f0 - 12.0f0)
        @test forcing.nh4_deposition == 2.0f0 .* forcing.no3_deposition
    end
end
