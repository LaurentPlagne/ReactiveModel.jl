# Solver backends

A backend owns the solver state of a `LinearProblem`. HiGHS is the default;
Clp is also included. Optional extensions are loaded only after their solver
package is installed and imported.

## Included backends

| Backend | Variable domains | Matrix updates | Export | Interface |
| :--- | :--- | :--- | :--- | :--- |
| `HiGHSBackend()` | Continuous, integer, binary | Yes | Yes | C API from `HiGHS_jll` |
| `ClpBackend()` | Continuous | Yes | Yes | C API from `Clp_jll` |

The HiGHS backend does not load HiGHS.jl or MathOptInterface. Comparison
benchmarks install those packages separately. A persistent backend retains the
numerical model between solves; basis reuse and refactorization remain solver
decisions. Integer solves do not provide LP row duals.

```@example backend_choice
using ReactiveModel

@model BackendExample(data, point) begin
    @variable x in [0.0, 10.0]
    @constraint floor x >= point.minimum
    @objective Min 2.0 * x
end
shared = ConstData(BackendExample(nothing))
point = (minimum=3.0,)
highs_problem = LinearProblem(shared, point; backend=HiGHSBackend(time_limit=10.0))
clp_problem = LinearProblem(shared, point; backend=ClpBackend())
highs_value = objective_value(highs_problem, point)
clp_value = objective_value(clp_problem, point)
@assert isapprox(highs_value, 6.0)
@assert isapprox(clp_value, highs_value)
(highs_value, clp_value)
```

HiGHS keywords are native solver options, with values checked by the backend.
RM configures simplex, sequential LP operation and logging during initialization;
constructor options do not override every setting applied at that stage.
`ClpBackend()` currently exposes no solver-option keywords.

Choose a solver by its required capabilities and verified behavior on your
workload. Numerical consistency on a few models does not establish a general
reliability or performance ranking.

## Optional extensions

| Backend | Load first | Current documented scope |
| :--- | :--- | :--- |
| `TinyHiGHSBackend()` | `using TinyHiGHS` | Continuous LPs, numerical matrix updates and export |
| `ClpJLBackend()` | `using ClpJL` | Continuous LPs; matrix updates and export raise `ArgumentError` |
| `CPLEXBackend()` | `using CPLEX` | Initial LP/MILP import; incremental updates and export raise `ArgumentError`; small Community smoke models checked |
| `XpressBackend()` | `using Xpress` | LP/MILP import, numerical updates and export implemented; licensed integration unqualified |

Install each solver package in a consumer environment according to its own
installation and licensing instructions before loading it. No adjacent solver
checkout is required by the core package or documentation build. An unavailable
extension raises an error identifying the package to load. Presence of an
extension is not a claim of validation for every solver version or operation.

For ClpJL, a changed matrix can be loaded through `rebuild!` before solving;
incremental matrix updates are unsupported. HiGHS or Clp supports that operation
directly. A failed update invalidates results and forces reconstruction before
the next attempt, rather than silently retaining the old matrix.

CPLEX also requires `rebuild!` when bounds, right-hand sides, coefficients or
costs change. These operations fail explicitly until incremental support is
implemented. Its supported constructor options are `silent` and `time_limit`;
other keywords raise `ArgumentError`. CPLEX and Xpress clones preserve their
constructor options. A time limit, interruption or unbounded model is reported
separately from a proof of infeasibility; an unknown native status remains an
error. These contracts do not establish that a licensed installation works.

## Shared assembly and solver policies

All six adapters consume the same mathematical assembly: column costs and
bounds, normalized row bounds, and a Julia sparse CSC matrix. Native adapters
then convert its indices to zero-based `Cint`; pure Julia adapters use its
one-based indices. Character encodings for row senses and variable domains are
shared by the commercial adapters. Native calls and result extraction remain
specific to each solver.

Import policies are explicit because removing a small coefficient changes the
numerical problem. Filtering occurs before sparse duplicate aggregation.

| Backend | Coefficients at import | Objective sense |
| :--- | :--- | :--- |
| Clp | All entries retained, including explicit zeros | Passed to the solver |
| HiGHS | Entries with `abs(a) < 1e-9` protected during import, then restored | Passed to the solver |
| TinyHiGHS | Entries with `abs(a) <= SimplexOptions().small_matrix_value` removed | Passed to the solver |
| ClpJL | Entries with `abs(a) <= 1e-12` removed | Max costs negated in the adapter |
| CPLEX, Xpress | Entries with `abs(a) <= 1e-12` removed | Passed to the solver |

The shared helpers run during construction or reconstruction. Repeated solves
continue to use each adapter's incremental operations and persistent solver
state. Sharing source alone does not prove that every policy shares identical
compiled code or reduces total precompilation time.

## Validation scope

The backend assembly revision passed 2,257 assertions with HiGHS, Clp and the
optional TinyHiGHS integration, plus 17 JET reports. ClpJL was absent: its
capability guards and the shared assembler were checked, but its full solver
integration was not revalidated.

CPLEX.jl 1.2.0 with the installed Community 22.1.1 library passed 109 additional
native smoke assertions on macOS ARM64: small LP/MILP models, min/max, LP duals,
rebuild, clone options, unsupported update recovery, infeasibility and
unboundedness, plus repeated LP/MILP imports on the same native handles. This
does not qualify larger models or other installations.
The Xpress.jl 0.18.2 extension was evaluated using its official switches that
skip native initialization; no Xpress solve was executed. License-free tests
cover capability guards and status classification. An optional solver's
absence or a successful definition load is not a successful integration check.

The package's native test application passed `--trim=safe` with the included
backends. An additional TinyHiGHS executable failed the trim verifier
on six dynamic `Val` calls inside the installed TinyHiGHS solver. The same
application and dependency failed at the same sites before this assembly
revision. TinyHiGHS's ordinary Julia integration passes; this dependency version
is not qualified for `--trim=safe`. ClpJL and the commercial adapters were not
qualified as trimmed executables.

The [synthetic benchmark workflow](parallel-results.md) identifies the backends,
wrappers and solver settings to compare. Generate and retain results from the
source revision and environment you actually run.

## Implementing a backend

The backend interface is experimental. [`AbstractBackend`](@ref) documents
required operations and result accessors. Read `src/Backends/HiGHSBackend.jl` and
`src/Backends/ClpBackend.jl` for implementations used by the current suite.

A backend must provide current numerical result buffers, report solver status
accurately and implement every operation it advertises. Unsupported operations
should fail explicitly. Check construction, updates, infeasibility/recovery,
cloning, rebuilding, export and domain rejection before documenting support.
Result accessors use normalized rows and `Float64` values; absent LP duals must
be represented consistently with `backend_has_duals`.

The default accessors read conventional native-backend fields. Alternative
storage can override `ReactiveModel.backend_status`, `ReactiveModel.backend_objective_value`,
`ReactiveModel.backend_primal_value` and `ReactiveModel.backend_dual_value`.
The [design guide](guide/design.md) explains storage and configuration ownership.
HiGHS clones retain options configured through RM.
