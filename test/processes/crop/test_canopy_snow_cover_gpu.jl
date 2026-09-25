using NeuralCrop
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

@testset "CUDA LPJmL crop snow cover suppresses absorbed PAR" begin
    cells = 4096
    for cft in (cft1, cft3)
        crop = init_crop(Float32, cells, CuArray)
        pet = init_pet(Float32, cells, CuArray)
        state = test_model_state(crop; pet)
        crop.state.canopy.lai .= 2.0f0
        pet.par .= 20.0f0
        snow = CuArray(repeat(Float32[0, 0.1], cells ÷ 2))

        if cft === cft3
            apar_crop_maize!(cft, state, pet, snow)
        else
            apar_crop!(cft, state, pet, snow)
        end
        synchronize()

        fpar = Array(crop.auxiliary.canopy.fpar)
        apar = Array(crop.auxiliary.canopy.apar)
        @test all(>(0.0f0), @view(fpar[1:2:end]))
        @test all(>(0.0f0), @view(apar[1:2:end]))
        @test all(iszero, @view(fpar[2:2:end]))
        @test all(iszero, @view(apar[2:2:end]))
    end
end
