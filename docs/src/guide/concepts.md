# Data and dependencies

ReactiveModel separates a mathematical definition, study inputs and the point
used for one solve. Understanding this separation is necessary for correct
caching and numerical updates.

## Definition, data and point

| Object | Role | What may change |
| :--- | :--- | :--- |
| `ProblemDefinition` | Variable/constraint families, keys, domains, relations and objective | Create a new definition when families or keys change |
| `ConstData` | Definition and study payload, held by reference | Payload contents may change between rebuilds |
| Parameter point | Selectors and numerical state for one solve | Values may change; the concrete `isbits` type stays fixed |
| `LinearProblem` | Mutable numerical workspace and backend | Updated by `solve`, `objective_value` and `rebuild!` |

```@raw html
<div class="rm-diagram" tabindex="0" role="group" aria-label="Definition, inputs and owned numerical state">
```

![The fixed definition and study payload are held by reference in ConstData. The changing point feeds one LinearProblem's private numerical workspace and persistent backend. Mutated study data requires rebuilding every sharing problem; changed families or keys require a new definition.](../assets/mechanism-model-state.svg)

```@raw html
</div>
```

Study data can contain arrays, dictionaries and mutable structs. Their lengths
need not be encoded in types. The parameter point can be a named tuple such as
`(t=1, scenario=2)` or an immutable struct of numbers. Ordinary arrays are not
`isbits` and belong in study data. A fixed-size numerical state can use a tuple.

ReactiveModel stores the preceding point by value. It does not copy study arrays
or monitor their mutations. For one problem, all points must have the same
concrete type, including field names and field types.

## Rebuild after study-data mutation

Keep the payload contents constant between rebuilds. After a mutation, call
[`rebuild!`](@ref) before the next solve, even if the parameter point changes.
Otherwise, cached components may still describe the old data.

```@example study_contract
using ReactiveModel

@model StudyContractDispatch(data, point) begin
    @variable production[g in 1:2] in [0.0, data.capacity[g]]
    @constraint balance production[1] + production[2] == data.demand[point.t]
    @objective Min 10.0 * production[1] + 20.0 * production[2]
end

data = (capacity=[10.0, 10.0], demand=[5.0, 9.0])
problem = LinearProblem(ConstData(StudyContractDispatch(data), data), (t=1,))
initial_cost = objective_value(problem, (t=1,))
@assert isapprox(initial_cost, 50.0)

data.capacity[1] = 2.0
rebuild!(problem)
new_cost = objective_value(problem, (t=1,))
@assert isapprox(new_cost, 80.0)
(initial_cost, new_cost)
```

Rebuilding keeps the definition and family keys, invalidates earlier results,
and reconstructs numerical data. Buffers may be reallocated. Before a first
solve, numerical construction can remain pending until a point is supplied.

If several problems share the payload, stop their solves and result queries,
mutate the data, then rebuild every problem before resuming. A generation field
in the point does not replace this requirement. To change family sizes, keys,
ranges, filters or domains, create a new definition and problem.

## Inferred dependencies

Start application development with `LinearProblem(...; safe_mode=true)`.
This mode reevaluates all numerical formulas on every solve, independently of
the inferred or explicit conditions below, and applies actual differences.
`update_report(problem)` records changes that a condition would have missed and
reevaluations without numerical change. The constructor default is `false` for
existing code. See [business functions and safe updates](../tutorials/business-functions.md)
for the recommended validation workflow and its limits.

The DSL examines each numerical formula and the setup assignments it uses.
`demand = data.demand[point.t]` makes a RHS using `demand` depend on `t`, including
through chains of aliases and tuple destructuring.

Analysis separates variable bounds, matrix coefficients, normalized RHS values
and the objective (including its constant). A component with no parameter
dependency stays invariant between rebuilds. Passing the whole point into a
helper uses a conservative whole-point comparison.

Helpers must depend on their explicit inputs. The analysis does not inspect
function bodies, detect global state changes or follow external mutations.
They can contain arbitrary deterministic calculations of `(study_data, point)`;
they do not need to be simple table accesses.
The [functional API](../tutorials/first-model-functional.md) has no syntax
analysis; its default predicates reevaluate components at distinct points.

## Affine normalization

Both sides of a constraint may contain scalar constants and variable terms.
For example, `point.a * (x + point.shift) >= point.rhs` is stored as
`point.a * x >= point.rhs - point.a * point.shift`. The matrix depends on `a`;
the normalized RHS depends on `a`, `shift` and `rhs`.

```@example affine_constants
using ReactiveModel

@model AffineConstantExample(data, point) begin
    demand = data.demand[point.t]
    @variable x in [0.0, 10.0]
    @constraint balance x + 1.0 == demand
    @objective Min x + data.fixed_cost[point.t]
end

data = (demand=[3.0, 6.0], fixed_cost=[5.0, 7.0])
problem = LinearProblem(ConstData(AffineConstantExample(data), data), (t=1,))
first_value = objective_value(problem, (t=1,))
solution = solve(problem, (t=2,))
@assert isapprox(first_value, 7.0)
@assert isapprox(primal(solution, :x), 5.0)
@assert isapprox(objective_value(solution), 12.0)
@assert isapprox(activity(solution, :balance), 5.0)
(first_value, primal(solution, :x), objective_value(solution))
```

Objective constants are included in public objective values. `activity` reports
the normalized variable expression, excluding constants from the original
relation. Here the normalized balance is `x == demand - 1.0`, so its activity is
`5.0` at the second point. See [DSL grammar](dsl.md) for supported affine forms.

## Choosing a reactivity condition

Use inference first. If an explicit condition is needed, choose a field name,
a nonempty tuple of field names, or a predicate with two named arguments:

| Form | Meaning at a distinct point |
| :--- | :--- |
| `update_when :t` | Reevaluate when field `t` changes |
| `update_rhs_when (:t, :scenario)` | Reevaluate when either field changes |
| `update_when ((previous, current) -> rule(previous, current))` | Reevaluate when the predicate returns `true` |
| `update_when false` | Keep invariant between rebuilds |
| `update_when true` | Reevaluate at each distinct point |

Generated comparisons and the whole-point cache use `isequal`: two `NaN` values
compare equal; `0.0` differs from `-0.0`. A custom predicate chooses its own
comparison semantics. With `safe_mode=false`, the whole-point cache skips
predicates at an unchanged point, even for `update_when true`. Safe mode
reevaluates the formulas at an unchanged point too.

An explicit condition replaces inference for that component. It must cover all
dependencies, including affine constants after normalization. A missing dependency
can leave stale numbers in the solver with `safe_mode=false`. Safe mode
reevaluates the component, corrects detected differences and records the miss.

```@example reactivity_conditions
using ReactiveModel

@model ExplicitReactivity(data, point) begin
    @variable x in [0.0, point.limit] update_when :limit
    @constraint floor point.a * (x + point.shift) >= point.rhs update_lhs_when :a update_rhs_when (:a, :shift, :rhs)
    @objective Min point.cost * x + point.offset update_when ((previous, current) -> !isequal((previous.cost, previous.offset), (current.cost, current.offset)))
end

p = (limit=10.0, a=1.0, shift=1.0, rhs=3.0, cost=2.0, offset=7.0)
problem = LinearProblem(ConstData(ExplicitReactivity(nothing)), p)
@assert isapprox(objective_value(problem, p), 11.0)
q = (limit=9.0, a=2.0, shift=0.5, rhs=7.0, cost=4.0, offset=9.0)
solution = solve(problem, q)
@assert isapprox(primal(solution, :x), 3.0)
@assert isapprox(objective_value(solution), 21.0)
objective_value(solution)
```

`update_lhs_when` and `update_rhs_when` apply only to constraints. A general
`update_when` supplies both; a specific clause overrides it regardless of order.
Repeated clauses and malformed field tuples or predicates are rejected.

Compatibility forms include `when`, `changed_when`, `lhs_when`, `rhs_when`,
unquoted fields, and implicit `o`/`n` expressions. The older `changed(expr)`,
`δ(expr)` and `diff(expr)` compare inferred dependencies of `expr`, not its
computed result. Prefer the forms above in new code.

## Families and composition

Keys in one family must have a concrete, uniform type and be unique. Integers,
symbols and tuples are supported. Use stable equality and hashing for custom
keys. Typed empty key collections are allowed and contribute no rows or columns.

Definitions expose variable and constraint collections that can be combined
with `merge`. Duplicate family names are rejected. Cross-definition references
must be declared with `@external_variables` in the DSL and supplied by the final
definition. This composition changes the definition; it is not an operation for
adding families to an already-built solver model.
