# DSL grammar and limits

`@model` lowers a deliberately bounded algebraic language to Julia callbacks.
It inspects syntax and setup dependencies; it does not analyze arbitrary Julia
programs or inspect the bodies of helper functions.

Start with [the first tutorial](../tutorials/first-model.md) for a complete
introduction. This page is a syntax reference; the tables show schematic forms
rather than independent programs.

## Model structure

The following complete example includes a setup alias and tuple destructuring.
Run it after [installation](../installation.md):

```@example grammar_reference
using ReactiveModel
@model Name(data, point) begin
    indices = data.indices
    scale, offset = (point.scale, point.offset)
    @variable x[i in indices] in [0.0, data.capacity[i]]
    @constraint floor[i in indices] scale * (x[i] + offset) >= point.rhs
    @objective Min sum(x[i] for i in indices)
end

data = (indices=1:2, capacity=[10.0, 10.0])
point = (scale=2.0, offset=1.0, rhs=6.0)
problem = LinearProblem(ConstData(Name(data), data), point)
solution = solve(problem, point)
@assert isapprox(objective_value(solution), 4.0)
@assert all(isapprox(primal(solution, :x, i), 2.0) for i in data.indices)
objective_value(solution)
```

The short form `@model Name begin ... end` uses the names `data` and `param`.
`Name()` and `Name(data)` also use the omitted default argument names.
A signature accepts at most two distinct argument names, without type annotations.
The generated factory takes study data: `definition = Name(data)`.

A model body contains declarations and setup assignments. Setup can bind a name
or destructure a tuple. Type annotations on setup bindings are unsupported. Bindings must be unique, use earlier
bindings, and never shadow a model argument or decision variable. Setup aliases
must describe data or parameters, rather than algebraic expressions containing
decision variables. They are copied into the callbacks that need them, so helpers
must be pure functions of their explicit inputs, or use private model memoization
that is deterministic and transparent to the mathematical result. Helpers may be
called again: setup does not guarantee one call globally. Invalidate memoized
values when their mathematical inputs change, including study-data changes
before `rebuild!`. Concurrent problems need independent mutable data and caches.
Mutating arrays, redefining a
binding, loops, conditional declarations, function definitions, nested macros,
and standalone statements are unsupported inside the model body. Define helper
functions outside the model.

Variable names must be unique, as must constraint names. At most one objective
is allowed; an omitted objective uses zero minimization. Family names and
iterator names cannot shadow model arguments or setup bindings.
`@external_variables x y` (also `@known_variables`) marks names supplied by another
definition for composition; these names must actually exist in the merged
problem. It does not create variables.

## Families and bounds

| Construct | Accepted form |
| :--- | :--- |
| Scalar family | `x`, `balance` |
| Indexed family | `x[i in data.I]`, `x[i ∈ data.I]`, `x[i = data.I]` |
| Several indices | `x[i in data.I, j in 1:i]` |
| Filter | `x[i in data.I; i > 0]`, `x[i in data.I, j in 1:i; i > j && j > 0]` |
| Bounds | `x in [lower, upper]`, `x in (lower, upper)` |
| Bounds from a helper | `x in bounds(data, point)` returning a pair |
| Inequality bounds | `lower <= x <= upper`, `x >= lower`, `x <= upper` |
| Free variable | `@variable x` |
| Domain | `domain=:continuous`, `domain=:integer`, `domain=:binary` |

Each iterator binds one symbol, and names cannot repeat within a family. Bare
index spaces (`x[I]`), destructured iterator names, empty indexing (`x[]`), and
multiple semicolon filters are rejected. Combine filter conditions with `&&`.
Subsequent ranges may depend on earlier indices.

Family ranges, filters, and domains are structural: they depend only on study
data, never on the parameter point or decision variables. Bounds can depend on
the point, but never on decision variables. A new range, filter, or domain needs
a new definition. `rebuild!` reloads numbers for the existing keys.

Bounds are converted to `Float64` when the numerical model is built or updated.
Both `x in [lower, upper]` and `x in bounds_helper(data, point)` accept numeric
values such as `Int`, `Float32` or mixed pairs. A helper can return `(0, 10)`
without converting it itself. Infinite bounds retain their usual meaning;
binary bounds are intersected with `[0.0, 1.0]` after conversion. The same rule
applies to functional `variable_bounds` callbacks, whose public return types
are preserved.

## Algebraic expressions

| Construct | Meaning |
| :--- | :--- |
| `x`, `x[i]`, `x[i,j]` | Declared scalar or keyed variable |
| Scalar Julia expression | Data, parameter fields, indexing and pure helper calls |
| `a + b`, `a - b`, `-a` | Addition, subtraction, unary negation |
| `coefficient * expression` | Product with exactly one variable-dependent factor |
| `expression / coefficient` | Division by a variable-independent scalar |
| `condition ? a : b` | Either affine branch; condition independent of variables |
| `sum(expression for i in I)` | Sum of affine expressions |
| `sum(expression for i in I, j in J if predicate)` | Several iterator bindings and one filter |
| `sum(sum(expression for j in J) for i in I)` | Nested sums |

An affine sum containing variables handles an empty generator as zero. Its
constant terms are preserved, including through filters and branches. A purely
scalar Julia call retains ordinary Julia semantics, including its behavior for
empty collections.

A constraint contains exactly one relation, `lhs <= rhs`, `lhs == rhs`, or
`lhs >= rhs`; `≤` and `≥` are accepted. Either side may contain variables and
constants. A `begin ... end` around the relation is accepted, with exactly one
expression. The objective similarly contains exactly one affine expression,
with sense `Min` or `Max` (also `min`, `minimize`, `max`, `maximize`).

For example, `a * (x + shift) >= rhs` is normalized to
`a * x >= rhs - a * shift`. `activity` reports the variable part `a * x`.
Objective constants are included in `objective_value`.

Products of decision variables, variable denominators, powers and other
nonlinear functions of decision variables are rejected at expansion. So are
strict relations, ranged constraints, assignment used as a relation, broadcasting
with decision variables, `dot` with decision variables, and the flattened
multi-`for` form `sum(x[i,j] for i in I for j in J)`. Use a comma-separated
iterator list or nested `sum` calls. Arbitrary statements, macros and closures
inside algebraic expressions are unsupported; move scalar calculations into
pure helpers.

### Composing objective fragments

An objective-only model can refer to variable families declared elsewhere:

```@example objective_composition
using ReactiveModel

@model Decisions(data, point) begin
    @variable x in [0.0, 10.0]
end

@model EnergyCost(data, point) begin
    @external_variables x
    @objective Min point.price * x + data.fixed_cost
end

@model PenaltyCost(data, point) begin
    @external_variables x
    @objective Min point.penalty * x + point.startup
end

data = (fixed_cost=2.0,)
base = Decisions(data)
objective = compose_objectives(EnergyCost(data).objective, PenaltyCost(data).objective)
definition = ProblemDefinition(variables=base.variables, constraints=base.constraints,
                               objective=objective)
problem = LinearProblem(ConstData(definition, data), (price=3.0, penalty=1.0, startup=5.0))
solution = solve(problem, (price=3.0, penalty=1.0, startup=5.0))
@assert objective_value(solution) == 7.0
solution = solve(problem, (price=3.0, penalty=-4.0, startup=2.0))
@assert objective_value(solution) == -6.0
```

All fragments must have the same sense. Emission follows argument order; costs
on the same columns and constants add, including cancellations. Function-based
`Objective` values can be mixed with DSL fragments. Composition only combines
objectives: supply variables and constraints in the final definition, and keep
their names, keys and dimensions fixed.

The composite predicate is the short-circuit OR of its children. When any child
requests an update, all fragments emit their complete current contribution
again. Safe mode checks the resulting aggregate objective and produces one
objective diagnostic, rather than a diagnostic for each fragment.

After changing study data, invalidate private memoization and call `rebuild!`
on each affected problem. New keys require a new definition. To give a clone
independent data, use `similar(problem; data=independent_data)` with the same
data type and fixed keys. Concurrent problems must own independent mutable
caches as well; closures capturing shared mutable state remain shared.

## Reactivity and validation

See [Choosing a reactivity condition](concepts.md#Choosing-a-reactivity-condition)
for explicit fields and predicates. Inference follows the setup aliases used by
each formula and keeps matrix coefficients separate from the normalized RHS.
It cannot detect a helper reading global state. Explicit conditions replace
inference and must cover all dependencies of the numeric component.

`test/test_dsl_grammar.jl` checks the accepted forms numerically with HiGHS and
Clp, compares incremental updates with fresh problems, and checks rejection of
unsupported declarations. The common tests/JET/trimming workload also exercises
filtered dependent families, nested sums and names that might otherwise capture
compiler temporaries. This validation covers the stated grammar and supported
contracts, rather than arbitrary Julia metaprograms.
