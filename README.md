# NeuralCrop.jl

[![NeuralCrop](https://github.com/yunan-l/NeuralCrop.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yunan-l/NeuralCrop.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![DOI](https://img.shields.io/badge/DOI-10.48550%2FarXiv.2512.20177-blue)](https://doi.org/10.48550/arXiv.2512.20177)
[![License: EUPL-1.2](https://img.shields.io/badge/license-EUPL--1.2-blue)](LICENSE)

NeuralCrop.jl is a differentiable hybrid crop-modelling framework that combines
mechanistic crop and soil processes with trainable neural components. It can
be trained with observational data and supports simulations ranging from
individual sites to regional grids on CPUs and GPUs.

The repository contains the Julia implementation used in the NeuralCrop
manuscript. The `Default` configuration contains no neural-network component.
The `Hybrid` configuration introduces trainable processes. Phenology,
management, soil water, carbon, nitrogen, and energy processes remain explicit.

## Main capabilities

- Differentiable seasonal crop simulations using Enzyme.jl.
- CPU and GPU execution through KernelAbstractions.jl.
- Explicit crop phenology, management, soil hydrology, carbon and nitrogen
  cycling, and mass-balance diagnostics.
- Trainable process representations for hybrid simulation.
- Site-scale and gridded simulation workflows.

## Related paper

Lin, Y., Bathiany, S., Badri, M., Gelbrecht, M., Hess, P., Groenke, B.,
Heinke, J., Müller, C., and Boers, N. (2025).
[NeuralCrop: Combining physics and machine learning for improved crop yield projections](https://doi.org/10.48550/arXiv.2512.20177).
*arXiv preprint arXiv:2512.20177*.

```bibtex
@article{lin2025neuralcrop,
  title={NeuralCrop: Combining physics and machine learning for improved crop yield projections},
  author={Lin, Yunan and Bathiany, Sebastian and Badri, Maha and Gelbrecht, Maximilian and Hess, Philipp and Groenke, Brian and Heinke, Jens and M{\"u}ller, Christoph and Boers, Niklas},
  journal={arXiv preprint arXiv:2512.20177},
  year={2025}
}
```

## Installation

NeuralCrop.jl currently targets Julia 1.10. Clone the repository and instantiate
its project environment:

```bash
git clone https://github.com/yunan-l/NeuralCrop.jl.git
cd NeuralCrop.jl
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```

An NVIDIA GPU and a working
[CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) installation are required only
for GPU execution and GPU-specific tests.

## Examples

The `examples/` directory provides initial conditions, daily climate forcing,
and a wheat simulation notebook.

Regional entry points are available under `scripts/`. Large forcing data and
trained checkpoints are not stored directly in this repository.

## Testing

Run the CPU test suite with:

```bash
julia --project=. test/runtests.jl
```

Run the Enzyme differentiability tests separately with:

```bash
julia --project=test -e 'using Pkg; Pkg.develop(Pkg.PackageSpec(path=".")); Pkg.instantiate()'
julia --project=test test/runtests_ad.jl
```

GPU tests require a functional CUDA device:

```bash
julia --project=. test/runtests_gpu.jl
```

## Development status

This is research software under active development. Interfaces and
configurations may change while the manuscript is under review. For exact
reproducibility, use a tagged release or a cited commit rather than the moving
`main` branch.

## Contributing

Questions, issue reports, and contributions are welcome through the
[GitHub issue tracker](https://github.com/yunan-l/NeuralCrop.jl/issues).

## Acknowledgements

NeuralCrop.jl was developed with support from the
[Earth System Modeling group](https://www.asg.ed.tum.de/esm/home/) at the
Technical University of Munich and the
[FutureLab on Artificial Intelligence](https://www.pik-potsdam.de/en/institute/departments/complexity-science/research/artificial-intelligence)
at the Potsdam Institute for Climate Impact Research. This work received
funding from the China Scholarship Council (grant agreement 202303250017) and
the Horizon Europe ClimTip project (grant agreement 101137601).

## License

NeuralCrop.jl is released under the
[European Union Public Licence v1.2](https://eupl.eu/1.2/en).
