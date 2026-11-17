# Benchmark method and reproduction

The benchmark scripts use deterministic synthetic LPs with explicit graph and
coefficient choices. Read the [input specification](benchmark-inputs.md),
[comparison code](benchmarks-jump.md) and [serial/parallel protocol](parallel-results.md)
before collecting or interpreting a run.

## Check numerical consistency

From the repository root, in a POSIX shell:

```sh
julia --startup-file=no --project=benchmarks -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=benchmarks benchmarks/test_reproduce.jl
julia --startup-file=no --project=benchmarks benchmarks/test_frontends.jl small medium large
rm_consistency=$(mktemp -d)
julia --startup-file=no --project=benchmarks benchmarks/reproduce.jl "$rm_consistency" 5
```

The small serial comparison uses JuMP, MOI and POI on six points per preset.
It checks objectives and feasibility before timing. The larger parallel harness
also recomputes objectives from primal values, checks every assigned point and
records native library identity. Equivalent models can have different optimal
primal or dual solutions, so vector equality is not required.

## Choose a matched comparison

Keep the solver, options, point sequence and requested result queries comparable.
The included paths distinguish flat-array and named JuMP formulations, guarded
and unconditional updates, default options, and an explicit Clp cache reset.
Record adapter reloads as part of the implementation being measured. Equal
native solver binaries do not imply equal model-state reuse or entry points.

If the solver changes too, report a complete frontend/solver comparison rather
than attributing the difference to syntax. Record row/column counts, update
policies and numerical checks alongside timings.

## State timing boundaries

The [parallel command](parallel-results.md#Run-the-benchmark) excludes construction,
compilation, validation, warmup and pre-sample GC from repeated-solve timings.
It includes numerical updates, optimization, objective queries and, in parallel,
team dispatch, completion and reduction. Serial and parallel runs cover the
same total point set with an independent model per chunk.

Report sample counts, medians and observed ranges. A range is not a confidence
interval. Julia allocation probes exclude native memory; warmed kernels and
complete team traversals have different boundaries. The supplied RM timing
paths use `safe_mode=false`; safe checking needs its own measurement.

## Keep a reproducible record

Retain raw output with the source revision and local changes, manifest, Julia
and package versions, backend options, input dimensions and seeds, machine and
thread configuration, and timing exclusions. The scripts write local TOML
records; do not mix records from different source or environment versions.

Run campaigns sequentially on an otherwise available machine. Background load,
CPU frequency and cache state can affect observations. Record depot cache state
and validate installation and platform coverage separately.

Fixed-cut LP timings do not measure cut insertion, rebuilds or full SDDP
training. The [performance guide](why/performance.md) gives application adoption
criteria, and the [compilation protocol](compilation-results.md) separates
first-use costs from repeated execution.
