"""
    UpdateCheck

One component's cumulative observations in safe mode. `family` is a variable
or constraint family name, or `:objective`; `component` is `:bounds`,
`:coefficients`, `:rhs` or `:objective`. `checks` counts completed comparisons,
`changes` counts comparisons with different numerical values, `missed` counts
changes that the predicate would have skipped, and `unnecessary` counts positive
predicates whose reevaluation produced no numerical change.

An unnecessary reevaluation is an observation at the visited points, not proof
that a dependency can be removed. Counts are per complete component/family,
not per scalar coefficient. Initial construction is not a comparison.
"""
struct UpdateCheck
    family::Symbol
    component::Symbol
    checks::Int
    changes::Int
    missed::Int
    unnecessary::Int
end

struct _UpdateSafety
    enabled::Bool
    checks::Vector{UpdateCheck}
    positions::Dict{Tuple{Symbol,Symbol},Int}
end
function _UpdateSafety(enabled::Bool)
    return _UpdateSafety(enabled, UpdateCheck[], Dict{Tuple{Symbol,Symbol},Int}())
end

function _record_update_check!(
    safety::_UpdateSafety, family::Symbol, component::Symbol, predicted::Bool, changed::Bool
)
    safety.enabled || return nothing
    key=(family, component)
    i=get(safety.positions, key, 0)
    if i!=0
        previous = safety.checks[i]
        safety.checks[i] = UpdateCheck(
            family,
            component,
            previous.checks + 1,
            previous.changes + changed,
            previous.missed + (!predicted && changed),
            previous.unnecessary + (predicted && !changed),
        )
        return nothing
    end
    push!(
        safety.checks,
        UpdateCheck(
            family,
            component,
            1,
            Int(changed),
            Int(!predicted && changed),
            Int(predicted && !changed),
        ),
    )
    safety.positions[key]=length(safety.checks)
    return nothing
end

"""
    update_report(problem)

Return a copy of the cumulative [`UpdateCheck`](@ref) observations collected
with `LinearProblem(...; safe_mode=true)`. Safe mode corrects missed updates;
inspect `missed` to find incomplete conditions, and `unnecessary` to identify
candidate evaluation costs to investigate. No conditions are rewritten.

Reports cover completed numerical-component comparisons, including comparisons
before a solve found infeasible. They exclude initial construction and rebuilds.
An unchecked problem returns an empty report. Mutating the returned vector does
not affect the problem. See [`reset_update_report!`](@ref) to begin a new sample.
"""
update_report(problem) = copy(problem.update_safety.checks)

"""
    reset_update_report!(problem)

Clear safe-mode observations without changing numerical values, solver state or
the current solution. Future solves continue checking components in safe mode.
"""
function reset_update_report!(problem)
    empty!(problem.update_safety.checks)
    empty!(problem.update_safety.positions)
    return problem
end
