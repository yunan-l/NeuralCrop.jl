# Runtime benchmark

Run from the repository root:

```sh
julia --project=. benchmark/runtime.jl
```

The benchmark reports compile-warmed CPU throughput in cell-days per second,
allocated bytes, and the projected memory for the canonical active-cell
domain. Record these values before and after runtime changes; hardware-specific
numbers are not committed as universal pass/fail thresholds.
