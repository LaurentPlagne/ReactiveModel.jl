# API reference

Use this page to look up constructors and callbacks after working through the
[first model](tutorials/first-model.md) or the
[functional tutorial](tutorials/first-model-functional.md). Names can be qualified
with `import ReactiveModel as RM`.

Numerical callbacks receive the study payload and parameter value. The public
contract is [data and dependencies](guide/concepts.md); result lifetime and
nonoptimal outcomes are described in [the solve cycle](guide/solve-cycle.md).
The API is experimental. Backend implementation details may change.

## Study, points and model

```@docs
ConstData
Parameters
ProblemDefinition
LinearProblem
Workspace
SolveError
```

## Declarative DSL

```@docs
@model
@coeff
```

## Mathematical families

```@docs
AbstractVariable
VariableTag
VariableBlock
lower_bound
upper_bound
variable_bounds
AbstractConstraint
ConstraintTag
ConstraintBlock
lhs_terms!
rhs_value
AbstractObjective
ObjectiveTag
Objective
compose_objectives
objective_terms!
```

## Change predicates

The DSL infers predicates from syntax unless an explicit condition replaces
inference. Function-based families use their supplied predicates. A predicate
that suppresses reevaluation must cover every dependency of its component.


```@docs
has_changed
lhs_has_changed
rhs_has_changed
```

## Safe update checks

Pass `safe_mode=true` when constructing a `LinearProblem` to reevaluate all
formulas independently of change predicates, including at identical points.
Actual numerical differences are applied and observations can be inspected
after a control-loop run. See the [business-function tutorial](tutorials/business-functions.md).
For composed objectives, predicates are combined with a short-circuit OR and
any requested update reevaluates every fragment. Safe mode compares and reports
the single aggregate objective; it does not identify a faulty fragment predicate.

```@docs
UpdateCheck
update_report
reset_update_report!
```

## Solving and results

`Solution` is a reused live view, not a retained record of one point. Copy values
you need before another solve. Failed solves and rebuilds invalidate live results.
`activity` uses normalized variable terms; LP duals are unavailable for MILPs.


```@docs
solve
rebuild!
ReactiveModel.objective_value
termination_status
is_optimal
primal
dual
activity
ReactiveModel.Solution
```

## Introspection and coefficients

```@docs
set_coeff!
iteration_space
constraint_sense
variable_block_keys
variable_block_offset
constraint_block_keys
constraint_block_offset
ReactiveModel.full_type_string
```

## Backends

Custom storage can override the qualified result methods below. These public
extension points are not exported into the user's namespace. See the
[design guide](guide/design.md).

```@docs
AbstractBackend
HiGHSBackend
TinyHiGHSBackend
ClpBackend
ClpJLBackend
CPLEXBackend
XpressBackend
ReactiveModel.backend_status
ReactiveModel.backend_objective_value
ReactiveModel.backend_primal_value
ReactiveModel.backend_dual_value
```
