# Business functions and safe updates

In an application, bounds and coefficients often come from functions of the
form `f(study_data, study_parameters)`. They can contain branches, physical
correlations and other calculations. They do not need to be table lookups.
They must be deterministic from these two explicit inputs; study data retains
its read-only contract between rebuilds.

This example minimizes the operating cost of an input supplying a requested
output. Temperature and availability affect its capacity; temperature and a
boost setting affect conversion efficiency. For each fixed point the constraint
is linear in the decision variable, even though the parameter formulas are
nonlinear. The coefficients below are illustrative, not a calibrated physical model.

## Write ordinary Julia functions

The complete runnable source is `examples/business_functions.jl`. The following
excerpts are read from that file when building this page.

```@eval
using Markdown, ReactiveModel
source=read(joinpath(pkgdir(ReactiveModel),"examples","business_functions.jl"),String)
fragment=split(split(source,"function available_capacity";limit=2)[2],"\n@model BusinessDispatch";limit=2)[1]
Markdown.parse("```julia\nfunction available_capacity"*fragment*"\n```")
```

The DSL does not analyze helper bodies. Passing the complete point to a helper
gives a conservative whole-point dependency. This supports complex functions
without reproducing their logic in the driver.

## Define the numerical formulas once

```@eval
using Markdown, ReactiveModel
source=read(joinpath(pkgdir(ReactiveModel),"examples","business_functions.jl"),String)
fragment=split(split(source,"@model BusinessDispatch";limit=2)[2],"\n# A deliberately";limit=2)[1]
Markdown.parse("```julia\n@model BusinessDispatch"*fragment*"\n```")
```

There is no separate implementation of capacity, efficiency or demand in an
update routine. A control loop only supplies the next point and requests a solve.

## Start with safe mode

Use `LinearProblem(...; safe_mode=true)` while developing and validating the
application. It can also remain enabled in production when exhaustive formula
evaluation is appropriate. It reevaluates every numerical component on every
solve call, including identical points, independently of inferred or explicit
change conditions. It compares against the current numerical model and sends
only actual differences to the solver. The solver workspace is retained.

```@example business_functions
using ReactiveModel
include(joinpath(pkgdir(ReactiveModel),"examples","business_functions.jl"))
using .BusinessFunctionsExample

result=BusinessFunctionsExample.run_study()
@assert length(result.values)==8
result.values
```

The driver contains nested loops, a skipped candidate and a boost branch. Each
visited point is checked against an analytically computable solution:
`input = requested_output / conversion_efficiency`. No numerical setter or
update ordering appears in the driver.

```@eval
using Markdown, ReactiveModel
source=read(joinpath(pkgdir(ReactiveModel),"examples","business_functions.jl"),String)
fragment=split(split(source,"function run_study";limit=2)[2],"\nfunction demonstrate_missed_update";limit=2)[1]
Markdown.parse("```julia\nfunction run_study"*fragment*"\n```")
```

## Read the report before optimizing conditions

[`update_report`](@ref) returns independent snapshots of cumulative observations,
one per family/component. Initial construction is excluded.

```@example business_functions
report=update_report(result.problem)
@assert all(c.missed==0 for c in report)
[(family=c.family,component=c.component,checks=c.checks,changes=c.changes,
  missed=c.missed,unnecessary=c.unnecessary) for c in report]
```

- `missed`: numerical values changed although the condition returned false.
  Safe mode applied the change. Investigate the condition before using the
  predicate-based path.
- `unnecessary`: the condition returned true but reevaluation found no numerical
  change. This identifies an observed evaluation cost, not a useless solver
  update: RM already suppresses unchanged numerical setters.

Here capacity changes once in seven comparisons, so six reevaluations find
unchanged bounds. Inspection of `available_capacity` suggests that a condition
depending on temperature and availability could avoid those evaluations.
An unchanged value can also be caused by clipping or a branch. Observations at
visited points do not prove that a dependency can be removed everywhere.

## An incomplete condition is corrected

Suppose a capacity condition only mentions `availability`, overlooking the
temperature calculation inside `available_capacity`. The source contains this
deliberately incomplete model:

```@eval
using Markdown, ReactiveModel
source=read(joinpath(pkgdir(ReactiveModel),"examples","business_functions.jl"),String)
fragment=split(split(source,"@model IncompleteCapacity";limit=2)[2],"\nstudy_data()";limit=2)[1]
Markdown.parse("```julia\n@model IncompleteCapacity"*fragment*"\n```")
```

```@example business_functions
corrected=BusinessFunctionsExample.demonstrate_missed_update()
bounds_check=only(c for c in update_report(corrected) if c.component===:bounds)
@assert bounds_check.missed==1
@assert objective_value(corrected,(temperature=35.0,availability=1.0,
    load_factor=0.7,boost=false,price=1.0))==27.0
bounds_check
```

Temperature rises from 20 to 35, availability stays at 1, and capacity falls
from 30 to 27. Safe mode applies the correct bound despite the incomplete
condition. The report records the missed update. The mode does not rewrite
conditions or silently claim that they are now valid.

## Choose the execution policy

Keep safe mode if the extra formula evaluation is acceptable. After reviewing
conditions and testing representative points, a separate predicate-based
workspace can be created with `similar(result.problem; safe_mode=false)`.
It shares study data but owns its solver, results and report. Safe mode is
inherited by `similar` unless explicitly changed. Use
`reset_update_report!(problem)` to begin collecting a fresh set of observations.

The constructor default remains `safe_mode=false` for existing applications.
The recorded benchmark campaigns used this predicate-based path. Safe mode
does extra evaluation and reporting; those timings and zero-allocation
observations do not describe its cost.

Safe mode ensures that change conditions cannot suppress numerical reevaluation
or forwarding of detected changes. It does not check the correctness of the
business functions, family topology, concurrent mutations or native solver
implementation. Call `rebuild!` after study-data mutations and create a new
definition when family keys or dimensions change. See
[data and dependencies](../guide/concepts.md).

From the repository root:

```sh
julia --startup-file=no --project=. examples/business_functions.jl
```
