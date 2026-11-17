function _domain_bounds(domain::Symbol, lower, upper)
    # All callback forms cross this boundary during both assembly and updates.
    lower, upper = Float64(lower), Float64(upper)
    if domain === :binary
        lower = max(lower, 0.0)
        upper = min(upper, 1.0)
        lower > upper &&
            throw(ArgumentError("binary variable bounds exclude both 0 and 1: ($lower, $upper)"))
    end
    return lower, upper
end

"""Shared numerical arrays used while assembling families, without model/callback types."""
struct _AssemblyBuffers
    lower_values::Vector{Float64}
    upper_values::Vector{Float64}
    row_blocks::Vector{Tuple{Symbol,Any}}
    row_coefficients::Vector{Dict{Int,Float64}}
    row_rhs::Vector{Float64}
    objective_coefficients::Dict{Int,Float64}
    backend_variable_blocks::Vector{_BackendVariableBlock}
    variable_block_offsets::Vector{Int}
    backend_constraint_blocks::Vector{_BackendConstraintBlock}
    constraint_block_offsets::Vector{Int}
end
function _assembly_buffers(workspace)
    return _AssemblyBuffers(
        workspace.lower_values,
        workspace.upper_values,
        workspace.row_blocks,
        workspace.row_coefficients,
        workspace.row_rhs,
        workspace.objective_coefficients,
        workspace.backend_variable_blocks,
        workspace.variable_block_offsets,
        workspace.backend_constraint_blocks,
        workspace.constraint_block_offsets,
    )
end

struct _BuildVariables{D,P}
    buffers::_AssemblyBuffers
    data::D
    point::P
end
@noinline function (job::_BuildVariables)(block::VariableBlock)
    b=job.buffers
    lower_bounds=Float64[]
    upper_bounds=Float64[]
    for key in block.keys
        lower, upper=variable_bounds(block.variable, job.data, job.point.value, key)
        lower, upper=_domain_bounds(block.domain, lower, upper)
        push!(lower_bounds, lower)
        push!(upper_bounds, upper)
        push!(b.lower_values, lower)
        push!(b.upper_values, upper)
    end
    offset=length(b.lower_values)-length(lower_bounds)
    push!(
        b.backend_variable_blocks,
        _BackendVariableBlock(
            block.domain, offset, lower_bounds, upper_bounds, zeros(length(lower_bounds))
        ),
    )
    push!(b.variable_block_offsets, offset)
    return nothing
end

struct _BuildConstraints{D,P,E}
    buffers::_AssemblyBuffers
    data::D
    point::P
    emitter::E
end
function _backend_row(coefficients, rhs, sense)
    indices=sort!(collect(keys(coefficients)))
    return _BackendConstraint(indices, Float64[coefficients[i] for i in indices], rhs, sense)
end
@noinline function (job::_BuildConstraints)(block::ConstraintBlock)
    b, emit=job.buffers, job.emitter
    constraints=_BackendConstraint[]
    for key in block.keys
        empty!(emit.coefficients)
        lhs_terms!(block.constraint, emit, job.data, job.point.value, key)
        coefficients=copy(emit.coefficients)
        rhs=Float64(rhs_value(block.constraint, job.data, job.point.value, key))
        push!(b.row_coefficients, coefficients)
        push!(b.row_rhs, rhs)
        push!(b.row_blocks, (block.name, key))
        push!(constraints, _backend_row(coefficients, rhs, block.sense))
    end
    offset=length(b.row_rhs)-length(constraints)
    push!(
        b.backend_constraint_blocks,
        _BackendConstraintBlock(offset, constraints, zeros(length(constraints))),
    )
    push!(b.constraint_block_offsets, offset)
    return nothing
end

function _build_direct_data!(const_data::ConstData, workspace::Workspace, pp)
    buffers=_assembly_buffers(workspace)
    _foreach_family(
        empty!,
        (
            buffers.lower_values,
            buffers.upper_values,
            buffers.row_blocks,
            buffers.row_coefficients,
            buffers.row_rhs,
            buffers.objective_coefficients,
            buffers.backend_variable_blocks,
            buffers.variable_block_offsets,
            buffers.backend_constraint_blocks,
            buffers.constraint_block_offsets,
        ),
    )
    workspace.backend_objective=nothing
    data=const_data.data
    _foreach_family(_BuildVariables(buffers, data, pp), const_data.definition.variables.blocks)
    emit=_Emitter(workspace.variable_indices, Dict{Int,Float64}())
    _foreach_family(_BuildConstraints(buffers, data, pp, emit), const_data.definition.constraints.blocks)
    empty!(emit.coefficients)
    emit.constant=0.0
    _emit_objective!(const_data.definition.objective, emit, const_data, pp)
    workspace.objective_offset[]=emit.constant
    merge!(workspace.objective_coefficients, emit.coefficients)
    indices=sort!(collect(keys(workspace.objective_coefficients)))
    objective=_BackendObjective(
        _BackendObjectiveCoefficients(
            indices, Float64[workspace.objective_coefficients[i] for i in indices]
        ),
        const_data.definition.objective.sense,
    )
    workspace.backend_objective=objective
    backend_set_verbose_flag!(workspace.backend; verbose=false)
    initialize!(workspace.backend, workspace.backend_variable_blocks, workspace.backend_constraint_blocks, objective)
    workspace.built=true
    workspace.parameter_snapshot=pp
    return nothing
end
