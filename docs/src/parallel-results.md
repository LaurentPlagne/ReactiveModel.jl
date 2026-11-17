# Serial and parallel benchmark protocol

The supplied benchmark compares repeated solves on deterministic synthetic
cascade LPs. The [input specification](benchmark-inputs.md) describes the graphs,
seeds and coefficient formulas. Use this protocol to collect results for your
checkout, Julia version and machine.

## Implementations

The harness includes direct RM/Clp_jll and RM/HiGHS_jll paths, plus JuMP models
using Clp.jl and HiGHS.jl/MOI. Flat-array and named JuMP formulations have explicit
update policies, described on the [comparison code page](benchmarks-jump.md).
The benchmark dependencies are separate from ReactiveModel's runtime dependencies.

Clp variants distinguish ordinary updates, an explicit inner-cache reset before
setters, and default solver options. A reset can still lead to a complete native
reload. The reset variant depends on the adapter's cache layout; inspect that
layout when changing versions. A shared solver binary does not guarantee equal
state reuse or equal native entry points.

## Run the benchmark

From the repository root, in a POSIX shell:

```sh
julia --startup-file=no --project=benchmarks -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=benchmarks -t 2 benchmarks/test_parallel_comparison.jl
rm_parallel=$(mktemp -d)
julia --startup-file=no --project=benchmarks -t 8 benchmarks/parallel_comparison.jl "$rm_parallel" 11 8 32 9
```

The arguments select a new or empty output directory, repetitions, parallel
chunks, scenarios and stock points per scenario. This example uses 32 scenarios
and nine stocks, or 288 LPs per traversal. Adjust the chunk count to the available
ThreadTeams pool. Serial and parallel runs cover the same total point set,
with one persistent model per chunk and one concave region. Chunk ownership
keeps mutable solver state independent of the physical thread executing it.

Run experiments sequentially on an otherwise available machine. Record external
load and the threading configuration; these scripts do not control CPU frequency
or establish a quiet machine by themselves.

## Numerical checks

Before timing, each assigned point is checked for objective agreement
(`atol=1e-6`, `rtol=1e-8`), objective recomputation from primal values and scaled
variable/row feasibility (`1e-7`). The harness requires RM and JuMP/Clp to load
the same native Clp library and records its binary hash.

Different optimal primal or dual vectors can be valid. These checks establish
consistency for the generated problems, not a mathematical certificate for an
arbitrary model. Retain failed checks when diagnosing an environment change.

## Timing and memory boundaries

Construction, compilation, validation, warmup and pre-sample GC are outside the
repeated-solve timer. Numerical updates, optimization and objective queries are
inside it. Parallel traversals also include team dispatch, completion and
reduction. Configuration order rotates between repetitions.

The script writes raw samples and metadata to `results.toml`. Report medians
and observed ranges with the sample count; observed ranges are not confidence
intervals. Parallel wall time divided by total LP count measures throughput,
not individual solve latency.

Allocation probes separate the warmed point-loop kernel from the complete team
traversal. Julia heap counts exclude native solver memory. GC observations are
process-wide and do not identify the cause of a timing difference. These RM
paths use `safe_mode=false`; safe checking and application callbacks need their
own measurements.

## Record and interpret a run

Keep the raw output with the source revision, local changes, package manifest,
backend options, input parameters, machine and thread configuration. Compare
implementations on matching point sequences and numerical checks before drawing
conclusions about runtime or scaling. The [measurement guide](why/performance.md)
explains how to choose an application comparison.

This workload keeps cuts and family keys fixed. Rebuilds, cut insertion and the
cost of a complete dynamic-programming or SDDP algorithm require separate
experiments. See the [benchmark method](benchmarks.md) for consistency checks
and the [compilation protocol](compilation-results.md) for first-use costs.
