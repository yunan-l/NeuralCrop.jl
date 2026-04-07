# NeuralCrop

[![Build Status](https://github.com/yunan-l/NeuralCrop.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yunan-l/NeuralCrop.jl/actions/workflows/CI.yml?query=branch%3Amain)


NeuralCrop is based on the LPJmL ("Lund-Potsdam-Jena managed Land") model, which explicitly simulates carbon, water, energy, and nitrogen flows for both natural vegetation and agricultural crops at 0.5° × 0.5° (latitude × longitude) spatial and daily temporal resolution. LPJmL is a state-of-the-art process-based GGCM and contributes to the global model intercomparison networks AgMIP and ISIMIP. It has been comprehensively evaluated at the global scale. In NeuralCrop, neural networks are embedded to replace or augment key biological processes of LPJmL that are uncertain or simplified but directly observable, such as photosynthesis and soil moisture dynamics, or to emulate processes that are heuristic, such as carbon allocation. The considered process-based components within NeuralCrop are fully differentiable to enable a seamless integration with neural networks. 
