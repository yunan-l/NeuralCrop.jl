# Test layout

The test tree mirrors the corresponding source tree under `src/`.

- Run the default suite with `julia --project=. test/runtests.jl`.
- Run one process test directly, for example
  `julia --project=. test/processes/crop/test_lambda_solver_c4.jl`.
- Files ending in `_gpu.jl` require a functional NVIDIA GPU and are not
  included by the default CPU-compatible test entry point.
- Run the Enzyme differentiability suite with
  `julia --startup-file=no --project=test test/runtests_ad.jl` after
  `julia --startup-file=no --project=test -e 'using Pkg; Pkg.develop(Pkg.PackageSpec(path=".")); Pkg.instantiate()'`.
