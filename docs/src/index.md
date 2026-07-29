# NeuralCrop.jl

NeuralCrop.jl is a framework for differentiable, global gridded crop modeling that combines process-based crop and soil dynamics with trainable neural components. It is designed for hybrid Earth-system modeling workflows where physical structure and machine-learning flexibility are both required.

The package targets daily time-step simulations over many grid cells and supports CPU/GPU workflows through Julia's array and accelerator ecosystem.

## Installation

NeuralCrop.jl is not yet registered in the Julia General registry. You can install it directly from GitHub by typing `]` in the Julia REPL and running:

```julia
pkg> add https://github.com/yunan-l/NeuralCrop.jl.git
```

For development and examples, clone the repository and instantiate the project:

```bash
git clone https://github.com/yunan-l/NeuralCrop.jl.git
cd NeuralCrop.jl
julia --project=. -e "import Pkg; Pkg.instantiate()"
```

!!! compat "Julia 1.10 is recommended"
    NeuralCrop.jl currently targets Julia 1.10.x.

## Quick start

A typical NeuralCrop workflow has four steps:

1. Load climate, crop management, soil, and LPJmL-derived initial-state data.
2. Select a crop functional type and initialize model state containers on CPU or GPU.
3. Run a daily crop simulation with either process-based or hybrid components.
4. Save or analyze model outputs such as yield, GPP, LAI, soil water, and soil carbon pools.

The repository includes a small demonstration dataset and a Jupyter notebook in `examples/` for wheat simulations over 10 grid cells.

```julia
using NeuralCrop

# 1. Activate package environment (outside Julia):
# julia --project=. -e "import Pkg; Pkg.instantiate()"
#
# 2. Run the demonstration workflow:
# open examples/wheat_rainfed.ipynb
#
# 3. For scripted runs, use the exported pipeline:
# - load forcing and initial-condition data
# - initialize states with init_states!
# - run daily_crop_C3! or daily_crop_C4!
# - inspect output fields (yield, lai, gpp, swc, litc, ...)
```

## Table of contents

### Introduction

```@contents
Pages = [
    "introduction/basic_concepts.md",
    "introduction/software_architecture.md",
    "introduction/mathematical_formulation.md",
]
Depth = 2
```

### Running NeuralCrop

```@contents
Pages = [
    "running/installation.md",
    "running/input_data.md",
    "running/initialization.md",
    "running/daily_simulations.md",
    "running/training.md",
]
Depth = 2
```

### Extending NeuralCrop

```@contents
Pages = [
    "extending/core_interfaces.md",
    "extending/adding_process_modules.md",
    "extending/hybrid_components.md",
]
Depth = 2
```

### Models

```@contents
Pages = [
    "models/crop_model.md",
    "models/soil_model.md",
    "models/hybrid_model.md",
]
Depth = 2
```

### Processes

```@contents
Pages = [
    "processes/climate/overview.md",
    "processes/crop/overview.md",
    "processes/soil/overview.md",
    "processes/hybrid/overview.md",
    "processes/neural_network/overview.md",
    "processes/utilities/overview.md",
]
Depth = 2
```

### Other

```@contents
Pages = [
    "contributing.md",
    "api_index.md",
    "references.md",
]
Depth = 2
```
