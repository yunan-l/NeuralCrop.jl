function output_training!(output::Output,
                          photos::Photos,
                          crop::Crop,
                          soil::Soil,
                          μ::NamedTuple,
                          σ::NamedTuple
)

    output.gpp = vcat(output.gpp, reshape((photos.agd .- μ.gpp) ./ (σ.gpp .- μ.gpp), (1, :)))
    output.lambda = vcat(output.lambda, reshape((photos.lambda .- μ.lambda) ./ (σ.lambda .- μ.lambda), (1, :)))
    output.vmax = vcat(output.vmax, reshape((photos.vmax .- μ.vmax) ./ (σ.vmax .- μ.vmax), (1, :)))
    output.vegc = vcat(output.vegc, reshape((crop.vegc .- μ.vegc) ./ (σ.vegc .- μ.vegc) , (1, :)))
    output.carbon_sum = vcat(output.carbon_sum, reshape(crop.carbon_sum, (1, :)))
    output.biomass = vcat(output.biomass, reshape(crop.biomass, (1, :)))
    output.growing_mask = vcat(output.growing_mask, reshape(crop.isgrowing, (1, :)))

    output.litc = vcat(output.litc, reshape((soil.litc .- μ.litc) ./ (σ.litc .- μ.litc), (1, :)))
    output.fastc = vcat(output.fastc, reshape((soil.fastc .- μ.fastc) ./ (σ.fastc .- μ.fastc), (1, :)))
    output.slowc = vcat(output.slowc, reshape((soil.slowc .- μ.slowc) ./ (σ.slowc .- μ.slowc), (1, :)))

    output.swc = vcat(output.swc, reshape((soil.swc .- μ.swc) ./ (σ.swc .- μ.swc), (1, :)))

end

function output_finetune!(output::Output,
                          photos::Photos,
                          crop::Crop,
                          soil::Soil,
                          μ::NamedTuple,
                          σ::NamedTuple
)

    output.gpp = vcat(output.gpp, reshape((photos.agd .- μ.gpp) ./ (σ.gpp .- μ.gpp), (1, :)))
    output.reco = vcat(output.reco, reshape((crop.resp .+ Zygote.dropgrad(photos.rd) .+ soil.rh .- μ.reco) ./ (σ.reco .- μ.reco), (1, :)))
    output.growing_mask = vcat(output.growing_mask, reshape(crop.isgrowing, (1, :)))
    output.swc = vcat(output.swc, reshape((soil.swc .- μ.swc) ./ (σ.swc .- μ.swc), (1, :)))

end


function output_yield!(output::Output,
                       crop::Crop,
                       day
)
    if day == 365
        output.yield = vcat(output.yield, reshape(max.(crop.yield, 0.0f0), (1, :)))
    end

end
