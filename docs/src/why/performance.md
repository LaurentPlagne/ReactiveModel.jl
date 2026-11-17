# Measuring repeated solves

Performance is a property of a workload, implementation and environment.
ReactiveModel aims to reduce repeated model assembly and avoid unnecessary
numerical updates. The effect must be measured for the application in question.

## A decision criterion

The [serial/parallel protocol](../parallel-results.md) describes a comparison on
synthetic LPs with direct RM adapters and explicit JuMP implementations. The
[comparison code](../benchmarks-jump.md) identifies the formulation and update
policies; the [input specification](../benchmark-inputs.md) defines the generated
cases. Use these as a starting point for an experiment relevant to your workload.

To decide whether to adopt RM, measure your current approach first. Check whether
numerical updates and result handling are a material fraction of total time.
Then run a matched implementation on representative points, including setup,
rebuilds and required queries. Compare validated numerical outputs before timing.
Count how many solves are needed to amortize any extra setup time and weigh the
observed saving against adaptation effort and missing capabilities. Small
same-solver differences can be irrelevant when the application spends most of
its time elsewhere.

For SDDP, [the training example](../tutorials/sddp.md) distinguishes repeated
state/noise solves from cut updates requiring rebuilds. Fixed-cut LP timings
do not measure total training time.

## Separate the execution paths

| Measurement | What it includes |
| :--- | :--- |
| Process startup | Loading Julia, packages and any compilation reached by the program |
| Construction and first solve | Definition/workspace construction, numerical assembly and initial solve; state whether compilation is included |
| Update at a new point | Dependency checks, formula evaluation, numerical updates and reoptimization |
| Same-point cache | Predicate-based reuse of the previous outcome without another optimization |
| Result extraction | Scalar queries or collection materialization, stated separately |
| Parallel throughput | Total elapsed time divided by completed problem count, with worker count recorded |

A cache measurement does not estimate the cost of a changed problem. Parallel
time per LP is a throughput normalization; it is not the latency of an
individual solve. Result extraction can be part of the application cost even
when excluded from a solver-focused timing.

## What influences runtime

Changing only a few bounds may require different work from changing matrix
coefficients. A retained LP model lets a solver reuse state when appropriate,
but neither a fixed matrix nor a persistent backend guarantees a fixed number
of pivots or the absence of refactorization. Degeneracy, conditioning, presolve,
solver options and point order can affect the outcome and work performed.

The HiGHS backend uses a direct C interface. This design choice alone does not
establish a speed advantage. A timing of a native solver call includes work
inside the solver; it cannot be interpreted as a general cost of `ccall`.

## Memory measurements

The benchmark's RM paths use `safe_mode=false`, the predicate-based execution
path. Safe mode reevaluates every component and collects diagnostics; its cost
must be measured separately. Start with the [safe workflow](../tutorials/business-functions.md)
to validate an application before refining update conditions for performance.

Julia allocation measurements count Julia heap allocations. They exclude native
solver allocations. Callbacks, changed matrix rows and result collections can
allocate; measure them on the required application path after warmup.

Record allocations on the same path as the runtime measurement, after warmup,
and report the Julia and dependency versions. An allocation count or absence of
observed GC does not by itself explain a runtime difference or parallel scaling.
Causal explanations require additional controlled experiments or profiling.

## Comparing frontends and solvers

Hold the solver and options fixed when studying frontend behavior. When changing
the solver as well, describe the comparison as a complete frontend/solver
configuration. Record numerical row/column counts because equivalent formulations
can reach the solver in different forms.

Check objective agreement and feasibility first. Different optimal primal or
dual solutions can be valid for a degenerate LP. A speed comparison is not a
substitute for numerical validation, and a few synthetic examples do not rank
solver reliability across applications.

## Reproduce and report

The [benchmark protocol](../benchmarks.md) gives the short campaign command,
validation thresholds, timing exclusions and output format. It is an installation
and numerical consistency baseline. A larger performance study should state its
question, problem sizes, sampling method and repetition count, and retain raw
samples and source versions.

Keep each run's raw records with its source and environment. Conclusions apply
to the tested workload and configuration.
