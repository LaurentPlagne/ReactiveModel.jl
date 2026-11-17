# Results and the solve cycle

`solve(problem, point)` and `objective_value(problem, point)` use the same
numerical update path. This guide describes its lifecycle and the lifetime of
returned results.

The sequence below describes the default predicate-based path. With
`safe_mode=true`, steps 3–4 instead reevaluate every numerical formula on every
call, including an identical point; only actual changes are sent to the backend.
An unchanged numerical model can still reuse its solved outcome. See the
[safe-mode tutorial](../tutorials/business-functions.md) for reports and examples.

## From point to solver outcome

1. Check that the point has the concrete `isbits` type expected by the problem.
2. On the first solve, assemble numerical components and initialize the backend.
3. At an identical point (`isequal`), reuse the previous outcome. A previous
   nonoptimal outcome raises `SolveError` again.
4. At a distinct point, evaluate dependency conditions, then the affected formulas.
   Compare numerical values using `isequal` and send changes to the backend.
5. Optimize if required. If the numerical problem is unchanged, reuse its outcome.

This comparison uses exact Julia equality semantics, not a floating-point
closeness tolerance. Solver feasibility tolerances are a separate matter.
An objective-constant change is retained in public results even when no native
coefficient changes.

A retained LP model permits solver state reuse. The solver decides whether to
reuse a basis, refactorize or take other numerical steps. MILP solves have a
different algorithmic lifecycle; LP warm-start claims do not describe the whole
integer solve.

```@raw html
<div class="rm-diagram" tabindex="0" role="group" aria-label="Safe and selective solve paths">
```

![Initial construction assembles the model, initializes the backend and solves. Later safe calls evaluate all formulas; selective calls reuse the outcome at an identical point or use conditions at a distinct point. Actual numerical differences reach the backend; an unchanged numerical model reuses its prior outcome.](../assets/mechanism-solve-flow.svg)

```@raw html
</div>
```

Reevaluating a formula does not itself send a setter. The comparison of its
numerical values decides that. Reusing an outcome can raise the previous
`SolveError` again; it does not always return a solution.

Construction can also happen in `rebuild!` or while cloning a built problem.
Optimization waits for a solve call. An outcome is reusable only after the
numerical model has actually been optimized.

## A live solution and explicit snapshots

`solve` returns a reusable [`ReactiveModel.Solution`](@ref). Later successful
solves, including calls through `objective_value(problem, point)`, update that
same object. A variable holding it does not retain an earlier solution.

```@example result_lifetime
using ReactiveModel

@model LifetimeExample(data, point) begin
    @variable x in [0.0, 10.0]
    @constraint floor x >= point.rhs
    @objective Min x
end
problem = LinearProblem(ConstData(LifetimeExample(nothing)), (rhs=2.0,))
solution = solve(problem, (rhs=2.0,))
saved_value = objective_value(solution)
saved_x = primal(solution, :x)
saved_variables = solution.variables

objective_value(problem, (rhs=4.0,))
@assert isapprox(saved_value, 2.0)
@assert isapprox(saved_x, 2.0)
@assert isapprox(saved_variables.x[0], 2.0)
@assert isapprox(objective_value(solution), 4.0)
(saved_value, objective_value(solution))
```

For an indexed family, create a vector such as
`[primal(solution, :production, g) for g in keys]` to retain values.
Scalar `primal`, `dual` and `activity` queries do not materialize result
collections. Accessing `solution.variables` or `solution.constraints` builds
snapshot collections on demand. Their caches are cleared on subsequent updates.

## Nonoptimal outcomes and errors

`solve` and `objective_value` raise [`SolveError`](@ref) when the backend does not
report an optimal solution. Inspect `err.status`; do not treat every status as
infeasibility. An unbounded problem or a solver limit can also prevent optimality.
A prior live solution is invalid after a failed solve or rebuild, and queries
raise an error until a successful solve refreshes it. This includes the scalar
properties `solution.objective_value` and `solution.termination_status`, as well
as their function forms. Previously materialized snapshot collections retain
their own numerical values.

```@example result_lifetime
status = try
    solve(problem, (rhs=20.0,))
    :unexpected_success
catch err
    err isa SolveError || rethrow()
    err.status
end
@assert status == :infeasible

recovered = solve(problem, (rhs=3.0,))
@assert isapprox(primal(recovered, :x), 3.0)
(status, objective_value(recovered))
```

```@raw html
<div class="rm-diagram" tabindex="0" role="group" aria-label="Live result and retained snapshot timeline">
```

![At rhs 2 the live solution and saved values equal 2. At rhs 4 the same live solution reports 4 while saved values stay 2. The infeasible rhs 20 invalidates live queries. Recovery at rhs 3 refreshes the same live object to 3; the saved scalar and variable snapshot still equal 2.](../assets/mechanism-result-lifetime.svg)

```@raw html
</div>
```

```@example result_lifetime
@assert recovered === solution
@assert isapprox(saved_value, 2.0)
@assert isapprox(saved_variables.x[0], 2.0)
```

For a nonoptimal solver outcome, the numerical model is retained at the failed
point; subsequent points update from that state. Basis reuse is backend-specific.
An exception during an update or backend call can leave numerical state unknown;
that path invalidates results and forces reconstruction on the next attempt.
The same policy applies when a callback or backend initialization interrupts
`rebuild!`: correct the cause, then the next solve assembles the model again
from the current study data and requested point. A rejected incompatible
definition, detected before reconstruction begins, leaves the prior model
and solution usable.

| Event | Numerical model for the next call | Previous live result |
| :--- | :--- | :--- |
| Nonoptimal solver outcome | Retained at the failed point | Invalid until successful solve |
| Exception during assembly, update or backend call | Reconstructed before the next attempt | Invalid until successful solve |
| Incompatible definition rejected before assembly | Prior model retained | Prior result remains usable |

Already materialized snapshots keep their values in all three cases.

## Duals and activities

For continuous LPs, `dual` reads the backend's row dual for the normalized
constraint. Interpret its sign with the objective sense and equation orientation.
A degenerate LP can have multiple valid dual optima. Check a sensitivity
interpretation against the formulation rather than assuming dual equality across
solvers.

`activity` evaluates the normalized variable expression. It omits constants
moved into the RHS; see [affine normalization](concepts.md#Affine-normalization).
For HiGHS MILP results, LP dual buffers contain `NaN`: LP row duals are not
available for the integer optimum.

## Independent concurrent problems

A `LinearProblem` owns mutable state. Do not solve, rebuild or query it concurrently
with another operation on the same problem. Use `similar(problem)` to obtain a
fresh numerical workspace sharing the study data. `copy` has the same behavior;
`deepcopy` is rejected for stateful solver resources.

To give the clone independent application data, pass
`similar(problem; data=worker_data)`. The new data must have the same concrete
type, keys and dimensions. A built source produces a clone assembled once,
using the new data from its first callback; an unbuilt source produces an unbuilt
clone. Prepare `worker_data` with the copies your application needs before
cloning.

HiGHS, ClpJL, CPLEX and Xpress clone factories preserve their options configured
through RM. A different configurable backend must supply its own clone factory;
the generic fallback uses its zero-argument constructor. Direct native option
changes outside RM are not tracked. An explicit backend can replace that choice,
for example `similar(problem; backend=HiGHSBackend(time_limit=10.0))`.
Keep shared data read-only throughout concurrent work. The
[reservoir tutorial](../tutorials/bellman-recursion.md) uses ThreadTeams with
a separate problem per chunk, without relying on thread IDs.
