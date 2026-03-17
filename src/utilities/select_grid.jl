function get_lonlatidx(file_path::String,
                       year::Int) #using crop_yield.nc as input
    
    lonlatidx_set = []
    ds = NCDataset(file_path, "r")

    longitude = ds["longitude"]
    latitude = ds["latitude"]
    yield_data = ds["crop_yield"]

    for i in axes(longitude, 1)
        for j in axes(latitude, 1)
            one_cell = yield_data[i, j, year]
            if (!ismissing(one_cell) && one_cell != 0)
                lonlat_idx = [i, j]
                push!(lonlatidx_set, lonlat_idx)
            end
        end
    end
    close(ds)
    
    return lonlatidx_set
end


function get_lonlat(file_path::String,
                    year::Int) #using crop_yield.nc as input
    
    lonlat_set = []
    ds = NCDataset(file_path, "r")

    longitude = ds["longitude"]
    latitude = ds["latitude"]
    yield_data = ds["crop_yield"]

    for i in axes(longitude, 1)
        for j in axes(latitude, 1)
            one_cell = yield_data[i, j, year]
            if (!ismissing(one_cell) && one_cell != 0)
                lonlat = [longitude[i], latitude[j]]
                push!(lonlat_set, lonlat)
            end
        end
    end
    close(ds)
    
    return lonlat_set
end
