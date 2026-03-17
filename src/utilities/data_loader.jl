using Lux, CUDA, LuxCUDA
const gdev = gpu_device()

# Function to load one point data
function OnePoint_one_dimension(file_path::String, 
                                variable::String, 
                                lonlatidx)
    ds = NCDataset(file_path, "r")
    
    i, j = lonlatidx
    
    dataset = Float32.(ds[variable][i, j, :])
    close(ds)
    
    return dataset
end

function OnePoint_dimensions(file_path::String, 
                             variable::String, 
                             lonlatidx)
    ds = NCDataset(file_path, "r")
    
    i, j = lonlatidx
    
    dataset = Float32.(ds[variable][i, j, :, :])
    close(ds)
    
    return dataset
end


function get_growing_time_range(lonlatidx_set,
                                 sdate,
                                 hdate,
                                 year) 
                        
    """
    Args:
        lonlatidx_set: the index of longitude and latitude (720*280)
        sdate: sowing date
        hdate: harvesting date
        year: time range of year
    """

    time_range = Vector{UnitRange{Int}}()

    for i in axes(lonlatidx_set, 1)
    
        lonidx, latidx = lonlatidx_set[i]
        
        if sdate[lonidx, latidx, year] > hdate[lonidx, latidx, year]
            time_range_ = 365*(year-2)+sdate[lonidx, latidx, (year-1)] : 365*(year-1)+hdate[lonidx, latidx, year]
            push!(time_range, time_range_)
        else
            time_range_ = 365*(year-1)+sdate[lonidx, latidx, year] : 365*(year-1)+hdate[lonidx, latidx, year]
            push!(time_range, time_range_)
        end
            
    end

    return time_range

end


function get_training_time_range(lonlatidx_set,
                                 sdate,
                                 hdate,
                                 year) 
                        
    """
    Args:
        lonlatidx_set: the index of longitude and latitude (720*280)
        sdate: sowing date
        hdate: harvesting date
        year: time range of year
    """

    time_range = Vector{UnitRange{Int}}()

    for i in axes(lonlatidx_set, 1)
    
        lonidx, latidx = lonlatidx_set[i]
        
        if sdate[lonidx, latidx, year] > hdate[lonidx, latidx, year]
            time_range_ = 365*(year-2)+sdate[lonidx, latidx, (year-1)]-1 : 365*(year-1)+sdate[lonidx, latidx, year]-2
            push!(time_range, time_range_)
        else
            time_range_ = 365*(year-1)+sdate[lonidx, latidx, year]-1 : 365*year+sdate[lonidx, latidx, year+1]-2
            push!(time_range, time_range_)
        end
            
    end

    return time_range

end

function get_time_point(lonlatidx_set,
                        sdate,
                        hdate,
                        year) 
                        
    """
    Args:
        lonlatidx_set: the index of longitude and latitude (720*280)
        sdate: sowing date
        hdate: harvesting date
        year: time range of year
    """

    time_point = Vector{Vector{Int}}()

    for i in axes(lonlatidx_set, 1)
    
        lonidx, latidx = lonlatidx_set[i]
        
        if sdate[lonidx, latidx, year] > hdate[lonidx, latidx, year]
            time_point_ = [365*(year-2)+sdate[lonidx, latidx, (year-1)], 365*(year-1)+hdate[lonidx, latidx, year], 365*(year-1)+sdate[lonidx, latidx, year]-2]
            push!(time_point, time_point_)
        else
            time_point_ = [365*(year-1)+sdate[lonidx, latidx, year], 365*(year-1)+hdate[lonidx, latidx, year], 365*year+sdate[lonidx, latidx, year+1]-2]
            push!(time_point, time_point_)
        end
            
    end

    return time_point

end

function load_nc_file_one_dimension(file_path::String, 
                                    variable::String,
                                    timerange::UnitRange
)

    ds = NCDataset(file_path, "r")

    dataset = ds[variable][:, :, timerange]

    close(ds)

    return dataset
end

function load_nc_file_dimensions(file_path::String, 
                                 variable::String,
                                 timerange::UnitRange
)

    ds = NCDataset(file_path, "r")

    dataset = ds[variable][:, :, :, timerange]

    close(ds)

    return dataset
end


function ExtractDay(file_path::String, variable::String, day_index, save_path::String)
    
    ds = NCDataset(file_path, "r")
    extract_band = ds[variable][:, :, day_index]
    
    latitudes = ds["latitude"][:]
    longitudes = ds["longitude"][:]
    times = ds["time"][day_index]
    
    new_dataset = NCDataset(save_path, "c")
    
    defDim(new_dataset, "latitude", length(latitudes))
    defDim(new_dataset, "longitude", length(longitudes))
    defDim(new_dataset, "time", length(times)) 

    defVar(new_dataset, "latitude", latitudes, ("latitude",), attrib = OrderedDict(
                   "units" => "degrees_north",
                   "long_name" => "latitude",
                   "standard_name" => "latitude",
                   "axis" => "Y"))
    defVar(new_dataset, "longitude", longitudes, ("longitude",), attrib = OrderedDict(
                   "units" => "degrees_east",
                   "long_name" => "longitude",
                   "standard_name" => "longitude",
                   "axis" => "X"))
    defVar(new_dataset, "time", times, ("time",), attrib = OrderedDict(
                   "units" => "days since 2010-1-1 0:0:0",
                   "calendar" => ds["time"].attrib["calendar"]))
    
    extract_var = defVar(new_dataset, variable, extract_band, ("longitude", "latitude", "time"), attrib = OrderedDict(
               "units" => ds[variable].attrib["units"],
               "missing_value" => ds[variable].attrib["missing_value"],
               "_FillValue" => ds[variable].attrib["_FillValue"]))
    
    extract_var[:, :, :] = extract_band
    
    close(new_dataset)
    close(ds)
    
end

# change date into index
function TimeToIdx(year, month, day, time)
    target_date = DateTimeNoLeap(year, month, day)
    time_idx = findall(x -> x == target_date, time)[1]
    
    return time_idx
end