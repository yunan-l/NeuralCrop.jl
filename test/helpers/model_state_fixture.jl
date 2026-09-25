_test_backend(reference) = values -> begin
    destination = similar(reference, eltype(values), size(values))
    copyto!(destination, values)
    return destination
end

"""Build a canonical lifecycle state over existing crop and soil test arrays."""
function test_model_state(
    crop::NeuralCrop.Crop,
    soil::NeuralCrop.Soil;
    managed_land = nothing,
    pet = nothing,
    climbuf = nothing,
    weather = nothing,
    output = nothing,
)
    cells = length(crop.state.canopy.lai)
    T = eltype(crop.state.canopy.lai)
    reference = crop.state.canopy.lai
    backend = _test_backend(reference)
    managed_land = isnothing(managed_land) ? NeuralCrop.init_managed_land(T, cells, backend) : managed_land
    pet = isnothing(pet) ? NeuralCrop.init_pet(T, cells, backend) : pet
    climbuf = isnothing(climbuf) ? NeuralCrop.init_climbuf(T, cells, backend) : climbuf
    weather = isnothing(weather) ? NeuralCrop.init_weather(T, cells, backend) : weather
    output = isnothing(output) ? NeuralCrop.init_output(T, cells, backend) : output
    return NeuralCrop.model_state(climbuf, crop, pet, soil, managed_land, weather, output)
end

function test_model_state(crop::NeuralCrop.Crop; kwargs...)
    T = eltype(crop.state.canopy.lai)
    cells = length(crop.state.canopy.lai)
    backend = _test_backend(crop.state.canopy.lai)
    soil = NeuralCrop.init_soil(T, cells, T.(NeuralCrop.soilparams.soildepth), backend)
    return test_model_state(crop, soil; kwargs...)
end

function test_model_state(soil::NeuralCrop.Soil; kwargs...)
    T = eltype(soil.water.storage)
    cells = size(soil.water.storage, 2)
    backend = _test_backend(soil.water.storage)
    crop = NeuralCrop.init_crop(T, cells, backend)
    return test_model_state(crop, soil; kwargs...)
end
