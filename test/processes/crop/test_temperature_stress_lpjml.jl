using NeuralCrop
using Test

function lpjml_temperature_stress(cft, temperature::T, daylength::T) where {T}
    daylength < T(0.01) && return zero(T)
    maximum_temperature = cft.path == 1 ? T(45) : T(55)
    temperature > maximum_temperature && return zero(T)
    temperature >= T(cft.temp_co2.high) && return zero(T)
    k1 = T(2) * log(one(T) / T(0.99) - one(T)) /
         (T(cft.temp_co2.low) - T(cft.temp_photos.low))
    k2 = (T(cft.temp_co2.low) + T(cft.temp_photos.low)) / T(2)
    k3 = log(T(0.99) / T(0.01)) /
         (T(cft.temp_co2.high) - T(cft.temp_photos.high))
    low = one(T) / (one(T) + exp(k1 * (k2 - temperature)))
    high = one(T) - T(0.01) * exp(k3 * (temperature - T(cft.temp_photos.high)))
    return low * high
end

@testset "Temperature stress matches LPJmL" begin
    for T in (Float32, Float64), cft0 in (cft1, cft3)
        cft = convert_precision(T, cft0)
        temperatures = T[10, 15, 20, 25, 45, 56]
        cells = length(temperatures)
        crop = init_crop(T, cells, identity)
        pet = init_pet(T, cells, identity)
        state = test_model_state(crop; pet)
        pet.daylength .= T(12)
        temp_stress(cft, pet, state, temperatures)
        expected = lpjml_temperature_stress.(Ref(cft), temperatures, pet.daylength)
        @test crop.auxiliary.photosynthesis.temperature_stress ≈ expected rtol = eps(T) * T(16)
    end
end
