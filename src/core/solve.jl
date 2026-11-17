# Native backends share construction/update operations and the dispatched
# result protocol. Alternative backends may override result access without
# reproducing the included backends' storage fields.
function _build!(lp::LinearProblem, pp)
    lp.built && _invalidate_model!(lp)
    try
        lp.numerical_updates.state.current[] = _parameters(pp)
        _execute(lp.numerical_updates, _EXEC_BUILD)
    catch
        _invalidate_model!(lp)
        rethrow()
    end
    return nothing
end

function _materialize_variable_group(backend, group::_SolutionGroup)
    blocks = map(group.blocks) do block
        values = [backend_primal_value(backend, block.offset+i) for i in eachindex(block.keys)]
        return SolutionVariableBlock(block.name, block.keys, block.positions, block.offset, values)
    end
    return _SolutionGroup(blocks)
end

function _materialize_constraint_group(backend, row_coefficients, group::_SolutionGroup)
    blocks = map(group.blocks) do block
        rows = (block.offset+1):(block.offset+length(block.keys))
        duals = [backend_dual_value(backend, index) for index in rows]
        activities = [
            sum(
                coefficient * backend_primal_value(backend, column) for
                (column, coefficient) in row_coefficients[index];
                init=0.0,
            ) for index in rows
        ]
        return SolutionConstraintBlock(block.name, block.keys, block.positions, block.offset, duals, activities)
    end
    return _SolutionGroup(blocks)
end

function _materialize_variables(solution::Solution)
    layout = getfield(solution, :layout).variables
    backend = getfield(solution, :backend)
    groups = map(group->_materialize_variable_group(backend, group), getfield(layout, :groups))
    return SolutionCollection(groups, getfield(layout, :names), getfield(layout, :positions))
end

function _materialize_constraints(solution::Solution)
    layout = getfield(solution, :layout).constraints
    backend = getfield(solution, :backend)
    rows = getfield(solution, :row_coefficients)
    groups = map(group->_materialize_constraint_group(backend, rows, group), getfield(layout, :groups))
    return SolutionCollection(groups, getfield(layout, :names), getfield(layout, :positions))
end

function _solution(
    lp::LinearProblem{CD,P,B,Workspace{P,B,V,C,I}}, status::Symbol
) where {CD,P,B,V,C,I}
    generation_token = lp.workspace.solution_generation
    objective = _current_objective(lp)
    return Solution{B,V,C,I}(
        _result_layout(lp),
        lp.backend,
        lp.variable_indices,
        lp.row_coefficients,
        generation_token,
        generation_token[],
        objective,
        status,
        nothing,
        nothing,
    )
end

function _refresh_solution!(solution::Solution, objective::Float64, status::Symbol)
    setfield!(solution, :generation, getfield(solution, :generation_token)[])
    setfield!(solution, :objective_value, objective)
    setfield!(solution, :termination_status, status)
    setfield!(solution, :variables_cache, nothing)
    setfield!(solution, :constraints_cache, nothing)
    return nothing
end

function _optimize_backend!(lp::LinearProblem)
    return backend_solve!(
        lp.backend,
        lp.backend_variable_blocks,
        lp.variable_block_offsets,
        lp.backend_constraint_blocks,
        lp.constraint_block_offsets,
    )
end

_backend_status(lp::LinearProblem) = backend_status(lp.backend)

_checked_parameters(::LinearProblem{C,P}, parameters::Parameters{P}) where {C,P} = parameters
_checked_parameters(::LinearProblem{C,P}, parameters::P) where {C,P} = Parameters(parameters)
function _checked_parameters(::LinearProblem{C,P}, parameters) where {C,P}
    throw(ArgumentError("parameters must have type $P"))
end

function _current_objective(lp::LinearProblem)
    return backend_objective_value(lp.backend) + lp.objective_offset[]
end

# Failure policy:
#
# A solve that ends without an optimum (infeasible, unbounded, ...) raises
# `SolveError` but does NOT discard the numerical model: every update for the
# failed point was applied before the solve, so the backend holds exactly that
# point's LP, and the next point is diffed against it as usual. Rebuilding
# instead would require a full model load and cold solve for each infeasible
# probe, even though its complete numerical model is already available.
#
# Whether the *basis* left by the failed solve is reused is the backend's
# decision, because it is a numerical-fidelity question, not a bookkeeping one:
# on degenerate LPs a warm start from it may reach another optimal vertex, hence
# `HiGHSBackend` keeps it; a
# Clp backend restarts cold. Objectives are unaffected either way: the simplex
# certifies any optimum it returns.
#
# Only an exception raised while updating or solving — a backend error, not an
# infeasibility — leaves the model in an unknown state; that path rebuilds.
function _invalidate_solution!(lp::LinearProblem)
    lp.workspace.solution_generation[] += UInt(1)
    solution = lp.solution
    if solution !== nothing
        setfield!(solution, :variables_cache, nothing)
        setfield!(solution, :constraints_cache, nothing)
    end
    return nothing
end

# A partial assembly/update or native exception cannot leave a reusable model.
# Keep this policy shared by explicit rebuilds and the solve recovery path.
function _invalidate_model!(lp::LinearProblem)
    _invalidate_solution!(lp)
    lp.parameter_snapshot = nothing
    lp.workspace.last_status = :not_solved
    lp.built = false
    return nothing
end

function _failed!(lp::LinearProblem, parameters)
    status = _backend_status(lp)
    status === :optimal && (status = :other_error)
    lp.parameter_snapshot = parameters
    lp.workspace.last_status = status
    throw(SolveError(status, "reactive model solve failed"))
end

function _prepare_backend!(lp::LinearProblem{C,P,B}, parameters::Parameters{P}) where {C,P,B}
    if !lp.update_safety.enabled &&
        lp.built &&
        lp.workspace.last_status !== :not_solved &&
        lp.parameter_snapshot !== nothing &&
        isequal(lp.parameter_snapshot, parameters)
        status = lp.workspace.last_status
        status === :optimal || throw(SolveError(status, "reactive model solve failed"))
        return :optimal
    end
    feasible = try
        needs_solve = true
        if lp.built
            snap = lp.parameter_snapshot
            if snap !== nothing
                snap_p = snap::Parameters{P}
                if lp.update_safety.enabled || !isequal(snap_p, parameters)
                    changed = _update_backend!(lp, snap_p, parameters, lp.backend)
                    if !changed && lp.workspace.last_status !== :not_solved
                        # The same LP as at the previous solved point: reuse its outcome.
                        # A prebuilt model still needs its first optimizer call.
                        lp.parameter_snapshot = parameters
                        lp.workspace.last_status === :optimal && return :optimal
                        needs_solve = false
                    end
                end
            end
        else
            _build!(lp, parameters)
        end
        needs_solve && first(_optimize_backend!(lp))
    catch
        # Any exception leaves unknown state, including a callback's SolveError.
        # Ordinary non-optimal solver outcomes are handled below this catch.
        _invalidate_model!(lp)
        rethrow()
    end
    _invalidate_solution!(lp)
    (feasible && _backend_status(lp) === :optimal) || _failed!(lp, parameters)
    lp.parameter_snapshot = parameters
    lp.workspace.last_status = :optimal
    solution = lp.solution
    solution === nothing || _refresh_solution!(solution, _current_objective(lp), :optimal)
    return :optimal
end

function objective_value(lp::LinearProblem{C,P,B}, parameters) where {C,P,B}
    parameters = _checked_parameters(lp, parameters)
    _prepare_backend!(lp, parameters)
    return _current_objective(lp)
end

"""Solve at `parameters` and return the reusable, lazily materialized live `Solution`."""
function solve(lp::LinearProblem{C,P,B}, parameters) where {C,P,B}
    parameters = _checked_parameters(lp, parameters)
    solution = lp.solution
    !lp.update_safety.enabled &&
        solution !== nothing &&
        isequal(lp.parameter_snapshot, parameters) &&
        lp.workspace.last_status === :optimal &&
        return solution
    status = _prepare_backend!(lp, parameters)
    solution = lp.solution
    solution === nothing || return solution
    solution = _solution(lp, status)
    lp.solution = solution
    return solution
end

"""
    variable_block_offset(lp::LinearProblem, name::Symbol) -> Int

Return the 0-based column offset of the variable block `name`.
"""
function variable_block_offset(lp::LinearProblem, name::Symbol)::Int
    return _result_layout(lp).variables[name].offset
end

"""
    constraint_block_offset(lp::LinearProblem, name::Symbol) -> Int

Return the 0-based row offset of the constraint block `name`.
"""
function constraint_block_offset(lp::LinearProblem, name::Symbol)::Int
    return _result_layout(lp).constraints[name].offset
end

"""
    variable_block_keys(lp::LinearProblem, name::Symbol)

Return the collection of keys for the variable block `name`.
"""
function variable_block_keys(lp::LinearProblem, name::Symbol)
    return _result_layout(lp).variables.by_name[name].keys
end

"""
    constraint_block_keys(lp::LinearProblem, name::Symbol)

Return the collection of keys for the constraint block `name`.
"""
function constraint_block_keys(lp::LinearProblem, name::Symbol)
    return _result_layout(lp).constraints.by_name[name].keys
end
