# Modeling approaches

Several Julia tools support repeated optimization. Their abstractions address
different needs; a choice should consider the model classes, solver support,
parameter contract and maintenance requirements of the application.

## Related tools

[JuMP](https://jump.dev/JuMP.jl/stable/) provides an algebraic modeling language
and access to multiple solvers. A model can be retained and modified between
solves: rebuilding at every point is not required. Its manual describes
[constraint modification](https://jump.dev/JuMP.jl/stable/manual/constraints/#Modify-a-constant-term).

[MathOptInterface](https://jump.dev/MathOptInterface.jl/stable/) (MOI) provides a
common solver interface, including model modification operations. Using it
directly gives applications explicit control over functions, sets and indices.
JuMP uses MOI to communicate with solvers.

[ParametricOptInterface](https://jump.dev/ParametricOptInterface.jl/stable/)
(POI) manages parameters in JuMP and MOI models. Its documentation shows how
to change a parameter and solve again. Its supported expressions and handling
of bounds should be checked for the package version in use.

ReactiveModel describes a fixed family of linear programs using callbacks or
its `@model` DSL. Numerical formulas receive study data and an `isbits`
parameter point. The DSL infers dependencies; explicit conditions can replace
that inference. The package retains a solver workspace and sends numerical
updates through its backend implementations.

These descriptions concern interfaces, not measured performance rankings.
The linked project manuals are the primary references for their capabilities.

## Choosing a contract

| Approach | How the application describes changes |
| :--- | :--- |
| JuMP model modifications | The application changes bounds, coefficients or constraint data through the modeling API. |
| POI parameters | The application updates parameter values represented in the formulation. |
| Direct MOI or solver API | The application manages model indices and modification calls. |
| ReactiveModel | Formulas depend on study data and a point; the library evaluates affected components. |

All require a correct mapping from application inputs to the numerical problem.
ReactiveModel's explicit predicates can omit dependencies if written incorrectly.
Its [safe mode](../tutorials/business-functions.md) reevaluates all formulas,
applies detected changes despite such predicates and reports misses. In the
predicate-based path, inference does not inspect helper bodies or monitor
external mutable state. Both modes retain the study-data and family-key contract.
A direct solver interface also requires maintaining backend-specific code.

## What the ReactiveModel DSL expresses

The [first tutorial](../tutorials/first-model.md) is a complete formulation and
execution example. It uses ordinary-looking algebra inside a bounded DSL:
`@model` does not accept arbitrary Julia programs. See the
[grammar reference](../guide/dsl.md) for accepted expressions and rejected forms.

The [functional API](../tutorials/first-model-functional.md) supplies the same
numerical components through callbacks. This can help when constructing model
families from application configuration. In both styles, changing dimensions
or family keys requires a new definition and problem.

## Comparing implementations

Use the same mathematical problem and verify objectives and feasibility before
comparing timings. Record both the formulation and the numerical model sent to
the solver: alternative parameter representations can change row counts and
presolve behavior. Update policies should be named explicitly, such as
unconditional updates or updates guarded by dependency checks.

The repository benchmark implementations are specific implementations of these
approaches. They do not represent all possible uses or optimizations of each
library. See the [benchmark protocol](../benchmarks.md) for reproduction and
[performance guide](performance.md) for interpretation.
