"""
TrainingDataLoader(data, data_index, device)

Build training inputs from forcing/parameter datasets.
"""
function TrainingDataLoader(data::NamedTuple, 
                            data_index::Vector{Int},
                            device;
                            training = false,
                            training_by_yield = false
)

    @unpack TrainingData = data

    if training && training_by_yield
        TrainingData = (
            output = TrainingData.output[:, data_index],   
            output_n = TrainingData.output_n[:, data_index],   
            yield = TrainingData.yield[:, data_index],   
            μ = TrainingData.μ,
            σ = TrainingData.σ,
            gdhy_yield = TrainingData.gdhy_yield[:, data_index]
        ) |> device
    else training
        TrainingData = (
            output = TrainingData.output[:, data_index],   
            output_n = TrainingData.output_n[:, data_index],   
            yield = TrainingData.yield[:, data_index],   
            μ = TrainingData.μ,
            σ = TrainingData.σ
        ) |> device
    end

    return TrainingData
end