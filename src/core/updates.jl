"""Updater inputs independent of the complete mathematical definition type."""
struct _UpdateContext{B,S}
    backend::B
    scratch::S
    buffers::_AssemblyBuffers
    objective_offset::Base.RefValue{Float64}
end

"""
Compute sparse coefficient changes into reusable buffers. Missing entries mean
numerical zero; removed nonzero entries therefore produce explicit zero setters.
The dictionaries remain complete rows/objectives and are not replaced by deltas.
"""
function _coefficient_changes!(
    indices, values, previous::Dict{Int,Float64}, current::Dict{Int,Float64}
)
    empty!(indices)
    empty!(values)
    for (column, value) in current
        isequal(get(previous, column, 0.0), value) && continue
        push!(indices, column)
        push!(values, value)
    end
    for (column, value) in previous
        (haskey(current, column) || isequal(value, 0.0)) && continue
        push!(indices, column)
        push!(values, 0.0)
    end
    return nothing
end

"""Stage one changed bound pair in owned model buffers and reusable setter scratch."""
@noinline function _stage_bounds!(
    native, column, local_index, lower_values, upper_values, domain, lower, upper, scratch
)
    lower, upper=_domain_bounds(domain, lower, upper)
    isequal(lower, lower_values[column]) && isequal(upper, upper_values[column]) && return false
    native.lower_bounds[local_index]=lower
    native.upper_bounds[local_index]=upper
    lower_values[column]=lower
    upper_values[column]=upper
    push!(scratch.columns, column)
    push!(scratch.lowers, lower)
    push!(scratch.uppers, upper)
    return true
end

"""Forward staged bounds through scalar or block setters, independently of callback types."""
@noinline function _commit_bounds!(backend, native, scratch)
    isempty(scratch.columns) && return nothing
    if length(scratch.columns)<=3
        for i in eachindex(scratch.columns)
            backend_change_col_bounds!(
                backend, scratch.columns[i], scratch.lowers[i], scratch.uppers[i]
            )
        end
    else
        update_variable_block_bounds!(backend, native, native.offset)
    end
    return nothing
end

"""
Apply one sparse row difference. Native row packets own their delta arrays:
shared scratch cannot be retained because evaluating the next row overwrites it.
The complete current row remains in `previous` for queries and future diffs.
"""
@noinline function _apply_matrix_row!(backend, native, row, previous, current, scratch)
    indices, values=scratch.changed_coefficients, scratch.changed_coefficient_values
    _coefficient_changes!(indices, values, previous, current)
    changed=!isempty(indices)
    if changed
        native.variable_indices=copy(indices)
        native.variable_coefficients=copy(values)
        update_constraint_variable_coefficients!(backend, row, native)
    end
    empty!(previous)
    merge!(previous, current)
    return changed
end

"""Forward only staged RHS rows, preserving the backend's batching policy."""
@noinline function _commit_rhs!(backend, native, first_row, rows)
    if length(rows)<=3
        for i in rows
            update_constraint_rhs!(backend, native.constraints[i], first_row+i-1)
        end
    else
        update_constraint_block_rhs!(backend, native, native.offset)
    end
    return nothing
end

struct VariableUpdater{P,K,V,CD,B,S} <: AbstractBlockUpdater
    backend::B
    data::CD
    block::VariableBlock{K,V}
    backend_blocks::Vector{_BackendVariableBlock}
    block_index::Int
    columns::UnitRange{Int}
    lower_values::Vector{Float64}
    upper_values::Vector{Float64}
    scratch::S
    updates::_UpdateState{P}
end

@noinline function _run_variable_update!(u::VariableUpdater{P}) where {P}
    block = u.block
    backend_block = u.backend_blocks[u.block_index]
    columns = u.columns
    updates = u.updates

    previous = updates.previous[]
    current = updates.current[]
    predicted = has_changed(block.variable, previous.value, current.value)
    (updates.safety.enabled || predicted) || return nothing
    numerical_change = false
    scratch = u.scratch
    empty!(scratch.columns)
    empty!(scratch.lowers)
    empty!(scratch.uppers)
    for (i, (key, column)) in enumerate(zip(block.keys, columns))
        lower, upper = variable_bounds(block.variable, u.data, current.value, key)
        numerical_change |= _stage_bounds!(
            backend_block,
            column,
            i,
            u.lower_values,
            u.upper_values,
            block.domain,
            lower,
            upper,
            scratch,
        )
    end
    updates.changed[] |= numerical_change
    _record_update_check!(updates.safety, block.name, :bounds, predicted, numerical_change)
    _commit_bounds!(u.backend, backend_block, scratch)
    return nothing
end

@inline (u::VariableUpdater)() = _run_variable_update!(u)

function _direct_variable_update(
    data, context, block::VariableBlock{K,V}, index::Int, first_column::Int, state::_UpdateState{P}
) where {K,V,P}
    columns = first_column:(first_column + length(block.keys) - 1)
    return VariableUpdater{P,K,V,typeof(data),typeof(context.backend),typeof(context.scratch)}(
        context.backend,
        data,
        block,
        context.buffers.backend_variable_blocks,
        index,
        columns,
        context.buffers.lower_values,
        context.buffers.upper_values,
        context.scratch,
        state,
    )
end

struct ConstraintUpdater{P,K,C,D,B,S} <: AbstractBlockUpdater
    backend::B
    data::D
    block::ConstraintBlock{K,C}
    backend_blocks::Vector{_BackendConstraintBlock}
    block_index::Int
    first_row::Int
    row_coefficients::Vector{Dict{Int,Float64}}
    row_rhs::Vector{Float64}
    scratch::S
    updates::_UpdateState{P}
end

@noinline function _run_constraint_update!(u::ConstraintUpdater{P}) where {P}
    block = u.block
    backend_block = u.backend_blocks[u.block_index]
    first_row = u.first_row
    updates = u.updates

    previous = updates.previous[]
    current = updates.current[]
    lhs_changed = lhs_has_changed(block.constraint, previous.value, current.value)
    rhs_changed = rhs_has_changed(block.constraint, previous.value, current.value)
    safe = updates.safety.enabled
    (safe || lhs_changed || rhs_changed) || return nothing
    numerical_lhs_change = false
    numerical_rhs_change = false
    emit = u.scratch.emitter
    changed_rows = u.scratch.rows
    empty!(changed_rows)
    for (i, key) in enumerate(block.keys)
        row = first_row + i - 1
        if safe || lhs_changed
            empty!(emit.coefficients)
            lhs_terms!(block.constraint, emit, u.data, current.value, key)
            numerical_lhs_change |= _apply_matrix_row!(
                u.backend,
                backend_block.constraints[i],
                row,
                u.row_coefficients[row],
                emit.coefficients,
                u.scratch,
            )
        end

        if safe || rhs_changed
            rhs = Float64(rhs_value(block.constraint, u.data, current.value, key))
            if !isequal(rhs, u.row_rhs[row])
                backend_constraint = backend_block.constraints[i]
                backend_constraint.rhs = rhs
                u.row_rhs[row] = rhs
                push!(changed_rows, i)
                updates.changed[] = true
                numerical_rhs_change = true
            end
        end
    end
    updates.changed[] |= numerical_lhs_change || numerical_rhs_change
    _record_update_check!(
        updates.safety, block.name, :coefficients, lhs_changed, numerical_lhs_change
    )
    _record_update_check!(updates.safety, block.name, :rhs, rhs_changed, numerical_rhs_change)
    _commit_rhs!(u.backend, backend_block, first_row, changed_rows)
    return nothing
end

@inline (u::ConstraintUpdater)() = _run_constraint_update!(u)

function _direct_constraint_update(
    data, context, block::ConstraintBlock{K,C}, index::Int, first_row::Int, state::_UpdateState{P}
) where {K,C,P}
    return ConstraintUpdater{P,K,C,typeof(data),typeof(context.backend),typeof(context.scratch)}(
        context.backend,
        data,
        block,
        context.buffers.backend_constraint_blocks,
        index,
        first_row,
        context.buffers.row_coefficients,
        context.buffers.row_rhs,
        context.scratch,
        state,
    )
end

struct ObjectiveUpdater{P,O,D,B,S} <: AbstractBlockUpdater
    backend::B
    data::D
    objective::Objective{O}
    objective_coefficients::Dict{Int,Float64}
    objective_offset::Base.RefValue{Float64}
    scratch::S
    updates::_UpdateState{P}
end

@noinline function _run_objective_update!(u::ObjectiveUpdater{P}) where {P}
    objective = u.objective
    updates = u.updates

    previous = updates.previous[]
    current = updates.current[]
    predicted = has_changed(objective.objective, previous.value, current.value)
    (updates.safety.enabled || predicted) || return nothing
    numerical_change = false
    emit = u.scratch.emitter
    empty!(emit.coefficients)
    emit.constant = 0.0
    _emit_objective_fragment!(objective.objective, emit, u.data, current.value)
    if !isequal(u.objective_offset[], emit.constant)
        u.objective_offset[] = emit.constant
        updates.changed[] = true
        numerical_change = true
    end
    coefficients = emit.coefficients
    old_coefficients = u.objective_coefficients
    changed_indices = u.scratch.changed_columns
    changed_values = u.scratch.changed_costs
    _coefficient_changes!(changed_indices, changed_values, old_coefficients, coefficients)
    if !isempty(changed_indices)
        numerical_change = true
        updates.changed[]=true
        backend_update_objective_coefficients!(u.backend, u.scratch.objective_coefficients)
    end
    empty!(old_coefficients)
    merge!(old_coefficients, coefficients)
    _record_update_check!(updates.safety, :objective, :objective, predicted, numerical_change)
    return nothing
end

@inline (u::ObjectiveUpdater)() = _run_objective_update!(u)

function _direct_objective_update(
    data, context, objective::Objective{O}, state::_UpdateState{P}
) where {O,P}
    return ObjectiveUpdater{P,O,typeof(data),typeof(context.backend),typeof(context.scratch)}(
        context.backend,
        data,
        objective,
        context.buffers.objective_coefficients,
        context.objective_offset,
        context.scratch,
        state,
    )
end

function _direct_numerical_updates(
    const_data::ConstData, workspace::Workspace, parameters::Parameters{P};
    state::_UpdateState{P}=_UpdateState(Ref(parameters), Ref(parameters), Ref(false), workspace.update_safety),
) where {P}
    return _make_numerical_updates(const_data, workspace, state)
end

# Keep plan construction out of the static execution adapter and its callers.
# Resolve keyword defaults above, so initial construction and rebuilds share
# this positional factory. Numerical update kernels keep their call boundaries.
@noinline function _make_numerical_updates(
    const_data::ConstData, workspace::Workspace, state::_UpdateState,
)
    context=_UpdateContext(
        workspace.backend,
        workspace.update_scratch,
        _assembly_buffers(workspace),
        workspace.objective_offset,
    )
    data=const_data.data
    variable_index = Ref(0)
    variable_offset = Ref(1)
    variables = _family_plan(const_data.definition.variables.blocks) do block
        variable_index[] += 1
        first_column = variable_offset[]
        variable_offset[] += length(block.keys)
        return _direct_variable_update(data, context, block, variable_index[], first_column, state)
    end
    constraint_index = Ref(0)
    constraint_offset = Ref(1)
    constraints = _family_plan(const_data.definition.constraints.blocks) do block
        constraint_index[] += 1
        first_row = constraint_offset[]
        constraint_offset[] += length(block.keys)
        return _direct_constraint_update(data, context, block, constraint_index[], first_row, state)
    end
    objective = _direct_objective_update(data, context, const_data.definition.objective, state)
    return _NumericalUpdates(state, variables, constraints, objective)
end

function _update_backend!(
    lp::LinearProblem{C,P,B}, old::Parameters{P}, new::Parameters{P}, ::B
) where {C,P,B<:AbstractBackend}
    updates = lp.numerical_updates
    state = updates.state
    state.previous[] = old
    state.current[] = new
    state.changed[] = false
    _execute(updates, _EXEC_UPDATE)
    return state.changed[]
end
