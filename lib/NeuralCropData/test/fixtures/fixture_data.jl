module FixtureData

using NCDatasets

function _coordinates!(ds, longitude, latitude)
    defDim(ds, "longitude", length(longitude))
    defDim(ds, "latitude", length(latitude))
    lon = defVar(ds, "longitude", Float64, ("longitude",))
    lat = defVar(ds, "latitude", Float64, ("latitude",))
    lon[:] = longitude
    lat[:] = latitude
    lon.attrib["units"] = "degrees_east"
    lat.attrib["units"] = "degrees_north"
    return nothing
end

function _write_climate_file(path, variable_name, units, longitude, latitude, offset)
    NCDataset(path, "c") do ds
        _coordinates!(ds, longitude, latitude)
        defDim(ds, "time", 6)
        time = defVar(ds, "time", Int32, ("time",))
        time[:] = Int32[2000, 2000, 2000, 2001, 2001, 2001]
        variable = defVar(
            ds, variable_name, Float32, ("latitude", "time", "longitude");
            fillvalue = -9999.0f0,
        )
        values = Array{Float32}(undef, length(latitude), 6, length(longitude))
        for longitude_index in eachindex(longitude), time_index in 1:6,
                latitude_index in eachindex(latitude)
            cell_value = 10 * latitude_index + longitude_index
            values[latitude_index, time_index, longitude_index] =
                Float32(offset + time_index + cell_value)
        end
        variable[:, :, :] = values
        variable.attrib["units"] = units
    end
    return path
end

function write_fixture(directory)
    longitude = [-0.5, 0.5, 1.5]
    latitude = [10.0, 11.0]
    cellid = Matrix{Union{Missing, Int32}}([
        2 0
        missing 3
        1 missing
    ])

    grid_path = joinpath(directory, "grid.nc")
    NCDataset(grid_path, "c") do ds
        _coordinates!(ds, longitude, latitude)
        variable = defVar(ds, "cellid", Int32, ("longitude", "latitude"); fillvalue = Int32(-9999))
        variable[:, :] = cellid
    end

    soil_path = joinpath(directory, "soil.nc")
    NCDataset(soil_path, "c") do ds
        _coordinates!(ds, longitude, latitude)
        # Reversed spatial dimension order exercises semantic alignment.
        soilcode = defVar(ds, "soilcode", Int32, ("latitude", "longitude"); fillvalue = Int32(-9999))
        soilph = defVar(ds, "soilph", Float32, ("longitude", "latitude"); fillvalue = -1.0f32)
        codes_lon_lat = Matrix{Union{Missing, Int32}}([
            9 1
            missing 14
            6 missing
        ])
        ph_lon_lat = Matrix{Union{Missing, Float32}}([
            8 6
            missing 9
            7 missing
        ])
        soilcode[:, :] = permutedims(codes_lon_lat)
        soilph[:, :] = ph_lon_lat
    end

    management_path = joinpath(directory, "management.nc")
    NCDataset(management_path, "c") do ds
        _coordinates!(ds, longitude, latitude)
        defDim(ds, "cft", 4)
        defDim(ds, "time", 2)
        cft = defVar(ds, "cft", Int32, ("cft",))
        time = defVar(ds, "time", Int32, ("time",))
        cft[:] = Int32[1, 2, 3, 4]
        time[:] = Int32[2000, 2001]

        dimensions = ("time", "latitude", "cft", "longitude")
        landuse = defVar(ds, "landfrac", Float32, dimensions; fillvalue = -9999.0f0)
        residue = defVar(ds, "residuefrac", Float32, dimensions; fillvalue = -9999.0f0)
        sowing_date = defVar(ds, "sdate", Float32, dimensions; fillvalue = -9999.0f0)
        phu = defVar(ds, "phusum", Float32, dimensions; fillvalue = -9999.0f0)
        fertilizer = defVar(ds, "fertilizer", Float32, dimensions; fillvalue = -9999.0f0)
        manure = defVar(ds, "manure", Float32, dimensions; fillvalue = -9999.0f0)
        landuse.attrib["units"] = "1"
        residue.attrib["units"] = "1"

        land_values = zeros(Float32, 2, 2, 4, 3)
        # CFT id 20. Coordinates are assigned by canonical cell id.
        land_values[1, 2, 2, 1] = 0.0f0 # cell 0
        land_values[1, 1, 2, 3] = 0.2f0 # cell 1
        land_values[1, 1, 2, 1] = 0.0f0 # cell 2
        land_values[1, 2, 2, 2] = 0.3f0 # cell 3
        land_values[2, 2, 2, 1] = 0.4f0
        land_values[2, 1, 2, 3] = 0.0f0
        land_values[2, 1, 2, 1] = 0.0f0
        land_values[2, 2, 2, 2] = 0.1f0
        land_values[:, :, 4, :] .= land_values[:, :, 2, :] .+ 0.1f0
        landuse[:, :, :, :] = land_values
        residue[:, :, :, :] = clamp.(land_values .+ 0.25f0, 0, 1)
        sowing_date[:, :, :, :] = ifelse.(land_values .> 0, 100.0f0, 0.0f0)
        phu[:, :, :, :] = ifelse.(land_values .> 0, 1200.0f0, 0.0f0)
        fertilizer[:, :, :, :] = 20.0f0 .* land_values
        manure[:, :, :, :] = 5.0f0 .* land_values
    end

    temp_path = _write_climate_file(
        joinpath(directory, "temp.nc"), "temp", "degC", longitude, latitude, -5,
    )
    prec_path = _write_climate_file(
        joinpath(directory, "prec.nc"), "prec", "mm/day", longitude, latitude, 0,
    )
    lwnet_path = _write_climate_file(
        joinpath(directory, "lwnet.nc"), "lwnet", "W/m2", longitude, latitude, -60,
    )
    swdown_path = _write_climate_file(
        joinpath(directory, "swdown.nc"), "swdown", "W/m2", longitude, latitude, 100,
    )
    co2_path = joinpath(directory, "co2.txt")
    open(co2_path, "w") do io
        write(io, "# year ppm\n2000 369.5\n2001 371.0\n")
    end

    catalog_path = joinpath(directory, "catalog.toml")
    open(catalog_path, "w") do io
        write(io, """
[cfts]
ids = [10, 20]
names = ["crop_a", "crop_b"]

[datasets.grid]
path = "grid.nc"
variable = "cellid"

[datasets.soilcode]
path = "soil.nc"
variable = "soilcode"

[datasets.soilph]
path = "soil.nc"
variable = "soilph"

[datasets.landuse]
path = "management.nc"
variable = "landfrac"
units = "1"
rainfed_bands = [1, 2]
irrigated_bands = [3, 4]

[datasets.residue_fraction]
path = "management.nc"
variable = "residuefrac"
units = "1"
rainfed_bands = [1, 2]
irrigated_bands = [3, 4]

[datasets.sowing_date]
path = "management.nc"
variable = "sdate"
rainfed_bands = [1, 2]
irrigated_bands = [3, 4]

[datasets.phu]
path = "management.nc"
variable = "phusum"
rainfed_bands = [1, 2]
irrigated_bands = [3, 4]

[datasets.fertilizer]
path = "management.nc"
variable = "fertilizer"
rainfed_bands = [1, 2]
irrigated_bands = [3, 4]

[datasets.manure]
path = "management.nc"
variable = "manure"
rainfed_bands = [1, 2]
irrigated_bands = [3, 4]

[datasets.temp]
path = "temp.nc"
variable = "temp"
units = "degC"

[datasets.prec]
path = "prec.nc"
variable = "prec"
units = "mm/day"

[datasets.lwnet]
path = "lwnet.nc"
variable = "lwnet"
units = "W/m2"

[datasets.swdown]
path = "swdown.nc"
variable = "swdown"
units = "W/m2"

[datasets.co2]
path = "co2.txt"
variable = "co2"
units = "ppm"
""")
    end
    return (;
        grid_path, soil_path, management_path, temp_path, prec_path, lwnet_path,
        swdown_path, co2_path, catalog_path,
    )
end

end
