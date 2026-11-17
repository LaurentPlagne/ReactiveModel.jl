# Frequently asked questions

## When is ReactiveModel relevant?

It targets repeated solves of linear programs with stable variable and constraint
families. Assess that contract alongside the solver capabilities, API stability
and application requirements. [Modeling approaches](why/ergonomics.md) describes
JuMP, MOI and POI with links to their own manuals. They also support model updates.

## Do I need change predicates?

For the DSL, start with inferred dependencies. It follows setup aliases and
parameter accesses but does not inspect helper bodies or observe global state.
Explicit conditions replace inference and must cover every dependency. The
functional API defaults to reevaluation at distinct points unless predicates
are supplied. See [data and dependencies](guide/concepts.md).

## What goes in the parameter point?

Small `isbits` values such as `(t=1, scenario=2)` or a tuple of numerical state.
Keep the same concrete type for all points of one problem. Dynamically sized
arrays belong in study data, whose contents stay constant between rebuilds.

## What if I change a study array?

Call `rebuild!` on every problem sharing that data before solving again, even
if the point also changes. RM does not detect external mutations. Stop concurrent
solves and queries during the change. New family keys or dimensions require
a new definition and problem. The [data guide](guide/concepts.md) includes an
executed mutation/rebuild example.

## Can a coefficient become zero or nonzero?

Numerical matrix coefficients can change, including zero/nonzero transitions,
with a backend that supports matrix updates. Variable and constraint keys stay
fixed. This differs from adding a variable or constraint family. Consult the
[backend capabilities](backends.md) and [DSL grammar](guide/dsl.md).

## How do I retain a result?

`solve` returns a reused live solution. Copy objective, primal or dual values
into your own scalars or arrays before the next solve. A failed solve or rebuild
invalidates earlier live results. The [solve-cycle guide](guide/solve-cycle.md)
shows retention and recovery with checked examples.

## What does `activity` include?

The normalized variable expression, excluding constants moved to the RHS.
For `x + 1 == demand`, activity is `x`, compared with `demand - 1`. Objective
constants are included in `objective_value`. See [normalization](guide/concepts.md#Affine-normalization).

## What happens if a solve is not optimal?

The call raises `SolveError`. Inspect its status: nonoptimal is not synonymous
with infeasible. The numerical model is retained after a nonoptimal solver
outcome; an exception during an update or backend call instead forces
reconstruction on the next attempt. See [errors](guide/solve-cycle.md).

## Can I solve concurrently?

Use one independent problem per task, with read-only shared data. `similar(problem)`
creates a new workspace; configure its backend explicitly if options matter.
`deepcopy` is rejected. The [reservoir tutorial](tutorials/bellman-recursion.md)
contains an executed ThreadTeams example with one workspace per chunk.

## Are integer and binary variables supported?

Yes, with HiGHS: `@variable x in [0, 1] domain=:binary` or
`@variable count in [0, 10] domain=:integer`. Domains stay fixed across solves.
Clp and TinyHiGHS reject these domains. LP row duals are unavailable for MILP
solutions. See [the domain example](guide/trimming.md#Integer-and-binary-variables).

## Is every warm solve allocation-free?

No. Allocation depends on callbacks, numerical updates, backend, Julia version
and result extraction. Native solver allocations are outside Julia's allocation
measurements. See [performance measurement](why/performance.md).

## What is `@coeff`?

A helper for callbacks implementing `lhs_terms!` or `objective_terms!`.
`@coeff emit production[g] = 1.0` emits a coefficient for a keyed variable.
It is not needed inside the algebraic `@model` DSL. The
[functional tutorial](tutorials/first-model-functional.md) uses explicit
`set_coeff!` calls, and the [API reference](api.md) describes extension points.

## How do I reproduce the examples and measurements?

[Installation](installation.md) provides the runnable example, package checks
and strict documentation build. [Benchmarks](benchmarks.md) gives a separate
numerical validation and measurement command. [Trimming](guide/trimming.md)
provides the compilation workload and its tested configuration.
