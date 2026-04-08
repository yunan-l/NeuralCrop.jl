<!-- Title -->
<h1 align="center">
NeuralCrop.jl
</h1>

<!-- description -->
<p align="center">
  <strong> 🧑‍🌾 💧 ☀️ 🌾 Fast and flexible Julia framework for hybrid crop modelling across scales. </strong>
</p>

<p align="center">
  <a href="https://yunan-l.github.io/NeuralCrop.jl/stable/">
    <img src="https://img.shields.io/badge/documentation-coming%20soon-orange" alt="Docs Status">
  </a>
  <a href="https://github.com/yunan-l/NeuralCrop.jl/actions">
    <img src="https://github.com/yunan-l/NeuralCrop.jl/workflows/CI/badge.svg" alt="Build Status">
  </a>
  <a href="https://doi.org/10.48550/arXiv.2512.20177">
    <img src="https://img.shields.io/badge/DOI-10.48550/arXiv.2512.20177-blue.svg" alt="DOI">
  </a>

</p>

NeuralCrop is a differentiable hybrid global gridded crop model (GGCM) that combines the strengths of the state-of-the-art GGCM [LPJmL](https://doi.org/10.5194/gmd-11-1343-2018) with machine learning approaches. By implementing process-based components in a differentiable form for seamless integration with machine learning methods, NeuralCrop enables end-to-end 'online training', with machine learning components optimized in tandem with the physical model dynamics. NeuralCrop is a flexible Julia framework supporting both purely process-based and hybrid simulations across CPUs and GPUs. More details is available in our preprint paper: [https://arxiv.org/abs/2512.20177](https://arxiv.org/abs/2512.20177)


## Installation

NeuralCrop is still under development to make it more user-friendly and not yet registered as a Julia package. To use it, you can still install the package from the repository via the package manager (type `]` in your REPL):
```
pkg> add https://github.com/yunan-l/NeuralCrop.jl.git
```

or clone the repository to your machine. 

Then, in the Julia REPL, activate the project and instantiate it to replicate our exact package versions:

```julia
pkg> activate(".")
pkg> instantiate()
```

This approach ensures you use the exact versions of all dependencies as specified in `Manifest.toml`, avoiding potential package version conflicts.

We recommend running NeuralCrop on Julia version 1.10.x.

## Example use

NeuralCrop does not provide the climate and management data required to drive the model, as these datasets originate from third-party sources. You can obtain the necessary input data from the [ISIMIP data repository](https://data.isimip.org/) (Inter-Sectoral Impact Model Intercomparison Project), and please cite the ISIMIP data appropriately when using it.

For a quick start, we provide a simplified demo in the examples/ directory, including 20-year forcing data (2000-2019) covering 10 grid cells. You can run the model in your Jupyter Notebook


## Citing

If you use NeuralCrop.jl in research or other activities 🏄, please mention NeuralCrop.jl and cite our paper:

> Lin, Yunan, et al. "NeuralCrop: Combining physics and machine learning for improved crop yield projections." arXiv preprint arXiv:2512.20177 (2025).

The bibtex entry for the paper is:

```bibtex
@article{lin2025neuralcrop,
  title={NeuralCrop: Combining physics and machine learning for improved crop yield projections},
  author={Lin, Yunan and Bathiany, Sebastian and Badri, Maha and Gelbrecht, Maximilian and Hess, Philipp and Groenke, Brian and Heinke, Jens and M{\"u}ller, Christoph and Boers, Niklas},
  journal={arXiv preprint arXiv:2512.20177},
  year={2025}
}
```