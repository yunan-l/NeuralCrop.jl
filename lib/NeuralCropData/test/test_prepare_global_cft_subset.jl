using Dates
using NCDatasets
using Test

include(joinpath(@__DIR__, "..", "scripts", "prepare_global_cft_subset.jl"))

@testset "Global CFT subset preparation" begin
    @test dimension_kind("pft") === :cft
    @test dimension_kind("cft") === :cft

    mktempdir() do directory
        management_path = joinpath(directory, "management.nc")
        NCDataset(management_path, "c") do dataset
            defDim(dataset, "longitude", 2)
            defDim(dataset, "band", 4)
            defDim(dataset, "latitude", 2)
            defDim(dataset, "time", 2)
            defVar(dataset, "longitude", Float64[0, 1], ("longitude",))
            defVar(dataset, "band", Int32[1, 2, 17, 18], ("band",))
            defVar(dataset, "latitude", Float64[1, 0], ("latitude",))
            defVar(dataset, "time", Int32[2000, 2001], ("time",))
            values = reshape(Float32.(1:32), 2, 4, 2, 2)
            defVar(
                dataset, "landfrac", values,
                ("longitude", "band", "latitude", "time"),
            )
        end
        management_output = joinpath(directory, "management_cft.nc")
        subset_netcdf(
            management_path, management_output, "landfrac";
            cft_index = 1,
            years = [2001],
            require_365_days = false,
            chunk_length = 1,
        )
        NCDataset(management_output, "r") do dataset
            @test size(dataset["landfrac"]) == (2, 1, 2, 1)
            @test dataset["band"][:] == Int32[1]
            @test dataset["landfrac"][:, 1, :, :] ==
                NCDataset(management_path, "r") do source
                    source["landfrac"][:, 1, :, 2:2]
                end
        end

        management_24_output = joinpath(directory, "management_24.nc")
        subset_netcdf(
            management_path, management_24_output, "landfrac";
            cft_indices = [1, 2, 1, 2], years = [2001],
            require_365_days = false, chunk_length = 1,
        )
        NCDataset(management_24_output, "r") do dataset
            @test size(dataset["landfrac"]) == (2, 4, 2, 1)
            @test dataset["band"][:] == Int32[1, 2, 3, 4]
            @test dataset.attrib["neuralcrop_source_cft_indices"] == "1,2,1,2"
            @test dataset["landfrac"][:, 1, :, :] == dataset["landfrac"][:, 3, :, :]
            @test dataset["landfrac"][:, 2, :, :] == dataset["landfrac"][:, 4, :, :]
        end

        pft_path = joinpath(directory, "fertilizer_pft.nc")
        NCDataset(pft_path, "c") do dataset
            defDim(dataset, "pft", 4)
            defDim(dataset, "time", 2)
            defVar(dataset, "pft", Int32.(1:4), ("pft",))
            defVar(dataset, "time", Int32[2000, 2001], ("time",))
            defVar(
                dataset, "fertilizer", reshape(Float32.(1:8), 4, 2),
                ("pft", "time"),
            )
        end
        pft_output = joinpath(directory, "fertilizer_2cfts.nc")
        subset_netcdf(
            pft_path, pft_output, "fertilizer";
            cft_indices = [2, 4], years = [2000, 2001],
            require_365_days = false, chunk_length = 1,
        )
        NCDataset(pft_output, "r") do dataset
            @test size(dataset["fertilizer"]) == (2, 2)
            @test dataset["pft"][:] == Int32[1, 2]
            @test dataset.attrib["neuralcrop_source_cft_indices"] == "2,4"
            @test dataset["fertilizer"][:, :] == NCDataset(pft_path, "r") do source
                source["fertilizer"][[2, 4], :]
            end
        end

        @test isnothing(management_years(management_path, "landfrac", nothing))
        management_all_output = joinpath(directory, "management_all.nc")
        subset_netcdf(
            management_path, management_all_output, "landfrac";
            cft_indices = [1, 2, 1, 2],
            years = management_years(management_path, "landfrac", nothing),
            require_365_days = false, chunk_length = 1,
        )
        NCDataset(management_all_output, "r") do dataset
            @test size(dataset["landfrac"]) == (2, 4, 2, 2)
            @test dataset["time"][:] == Int32[2000, 2001]
        end

        sdate_path = joinpath(directory, "sdate.nc")
        NCDataset(sdate_path, "c") do dataset
            defDim(dataset, "longitude", 2); defDim(dataset, "band", 24)
            defDim(dataset, "latitude", 2); defDim(dataset, "time", 1)
            defVar(dataset, "longitude", Float64[0, 1], ("longitude",))
            defVar(dataset, "band", Int32.(1:24), ("band",))
            defVar(dataset, "latitude", Float64[1, 0], ("latitude",))
            defVar(dataset, "time", Int32[0], ("time",))
            defVar(dataset, "sdate", reshape(Float32.(1:96), 2, 24, 2, 1),
                ("longitude", "band", "latitude", "time"))
        end
        @test has_time_dimension(sdate_path, "sdate")
        @test isnothing(management_years(sdate_path, "sdate", 2015))
        sdate_output = joinpath(directory, "sdate_24.nc")
        subset_netcdf(sdate_path, sdate_output, "sdate";
            cft_indices = collect(1:24), years = management_years(sdate_path, "sdate", 2015),
            require_365_days = false, chunk_length = 1)
        NCDataset(sdate_output, "r") do dataset
            @test size(dataset["sdate"]) == (2, 24, 2, 1)
        end

        management_only_config = joinpath(directory, "management_only.toml")
        open(management_only_config, "w") do output
            TOML.print(output, Dict(
                "subset" => Dict(
                    "output_directory" => joinpath(directory, "management_only"),
                    "management_years" => "all",
                    "chunk_length" => 1,
                ),
                "management" => Dict(
                    "phu" => Dict(
                        "input" => sdate_path,
                        "output" => "phu_24cfts.nc",
                        "variable" => "sdate",
                        "cft_indices" => collect(1:24),
                    ),
                ),
            ))
        end
        prepare_subset(management_only_config)
        @test isfile(joinpath(directory, "management_only", "phu_24cfts.nc"))

        climate_path = joinpath(directory, "climate.nc")
        time_values = Int32.(0:729)
        NCDataset(climate_path, "c") do dataset
            defDim(dataset, "time", length(time_values))
            defDim(dataset, "latitude", 1)
            defDim(dataset, "longitude", 2)
            defVar(
                dataset, "time", time_values, ("time",);
                attrib = Dict("units" => "days since 2001-01-01", "calendar" => "365_day"),
            )
            defVar(dataset, "latitude", Float64[0], ("latitude",))
            defVar(dataset, "longitude", Float64[0, 1], ("longitude",))
            values = reshape(
                Float32.(1:(length(time_values) * 2)), length(time_values), 1, 2,
            )
            defVar(dataset, "temp", values, ("time", "latitude", "longitude"))
        end
        NCDataset(climate_path, "r") do dataset
            _, indices = daily_indices_for_years(dataset, "temp", [2001, 2002], 2001)
            @test indices == collect(1:730)
        end
        climate_output = joinpath(directory, "climate_two_years.nc")
        subset_netcdf(
            climate_path, climate_output, "temp";
            years = [2001, 2002], daily_source_start_year = 2001, chunk_length = 31,
        )
        NCDataset(climate_output, "r") do dataset
            @test size(dataset["temp"]) == (730, 1, 2)
            @test calendar_year(first(dataset["time"][:])) == 2001
            @test calendar_year(last(dataset["time"][:])) == 2002
        end

        co2_path = joinpath(directory, "co2.txt")
        write(co2_path, "2000 368\n2001 370\n2002 372\n2003 374\n")
        co2_output = joinpath(directory, "co2_subset.txt")
        subset_co2(co2_path, co2_output, [2001, 2002])
        @test occursin("2001 370", read(co2_output, String))
        @test !occursin("2000 368", read(co2_output, String))
    end
end
