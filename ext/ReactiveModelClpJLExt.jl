module ReactiveModelClpJLExt

using ReactiveModel
import ReactiveModel: AbstractBackend, ClpJLBackend,
    _BackendVariableBlock, _BackendConstraintBlock, _BackendConstraint,
    _BackendObjective, _BackendObjectiveCoefficients, _backend_constraint_bounds, _backend_data,
    backend_has_duals, backend_set_verbose_flag!, backend_set_objective_sense!,
    backend_write_to_file, initialize!, backend_change_col_bounds!,
    backend_change_cols_bounds!, update_variable_block_bounds!,
    update_constraint_rhs!, update_constraint_block_rhs!,
    update_constraint_variable_coefficients!, backend_update_objective_coefficients!,
    backend_solve!

using ClpJL

function ReactiveModel.ClpJLBackend(; sparse_factorization::Union{Nothing, Bool}=nothing,
                                    options::ClpJL.SimplexOptions=ClpJL.SimplexOptions(presolve=false, scaling=false))
    ClpJLBackend(
        nothing, sparse_factorization, options, false,
        Float64[], Float64[], Float64[], Float64[], NaN, :not_solved,
    )
end

ReactiveModel._similar_backend(backend::ClpJLBackend) =
    ReactiveModel.ClpJLBackend(sparse_factorization=backend.sparse_factorization, options=backend.options)

backend_has_duals(::ClpJLBackend) = true

backend_set_verbose_flag!(::ClpJLBackend; verbose::Bool) = nothing

function backend_set_objective_sense!(backend::ClpJLBackend, sense::Symbol)
    model = backend.model::ClpJL.ClpModel
    old_dir = model.direction
    new_dir = sense === :max ? ClpJL.Maximize : ClpJL.Minimize
    if old_dir != new_dir
        model.direction = new_dir
        for i in eachindex(model.objective)
            model.objective[i] = -model.objective[i]
            model.workspace.objective_coefficients[i] = -model.workspace.objective_coefficients[i]
        end
        backend.objective_changed = true
    end
    nothing
end

function _crash_free_variables!(model::ClpJL.ClpModel)
    matrix = model.matrix
    basis = model.basis
    num_cols = model.num_cols
    num_rows = model.num_rows

    used_rows = falses(num_rows)
    for col in 1:num_cols
        if model.col_lower[col] <= -1e20 && model.col_upper[col] >= 1e20
            start_k = matrix.col_starts[col]
            end_k = matrix.col_starts[col + 1] - 1
            for k in start_k:end_k
                row = matrix.row_indices[k]
                if !used_rows[row] && abs(matrix.elements[k]) > 1e-4
                    slack_var = num_cols + row
                    basis.basic_variables[row] = col
                    basis.variable_statuses[col] = ClpJL.VarBasic
                    basis.variable_statuses[slack_var] = ClpJL.VarAtLowerBound
                    used_rows[row] = true
                    break
                end
            end
        end
    end
    nothing
end

function initialize!(backend::ClpJLBackend,
                     variable_blocks::Vector{_BackendVariableBlock},
                     constraint_blocks::Vector{_BackendConstraintBlock},
                     objective::_BackendObjective)
    all(block -> block.domain === :continuous, variable_blocks) ||
        throw(ArgumentError("ClpJLBackend supports continuous variables only"))

    data = _backend_data(variable_blocks, constraint_blocks, objective; drop_tolerance=1e-12)
    ncol, nrow = data.ncols, data.nrows
    direction = objective.op === :min ? ClpJL.Minimize : ClpJL.Maximize
    if direction == ClpJL.Maximize
        data.cost .*= -1.0
    end

    use_sparse = backend.sparse_factorization === nothing ? (nrow > 40) : backend.sparse_factorization
    factorization = use_sparse ?
        ClpJL.ClpSparseFactorization(nrow) : ClpJL.ClpDenseFactorization(nrow)

    model = ClpJL.ClpModel(factorization; options=backend.options)
    ClpJL.load_problem!(model, data.matrix, data.lower, data.upper, data.cost,
                        data.row_lower, data.row_upper; direction=direction)
    _crash_free_variables!(model)

    backend.model = model
    backend.objective_changed = false
    backend.col_value = zeros(ncol)
    backend.col_dual = zeros(ncol)
    backend.row_value = zeros(nrow)
    backend.row_dual = zeros(nrow)
    backend.objective_value = NaN
    backend.status = :not_solved
    nothing
end

function backend_change_col_bounds!(backend::ClpJLBackend, index::Int,
                                    lower::Float64, upper::Float64)
    ClpJL.set_column_bounds!(backend.model::ClpJL.ClpModel, index, lower, upper)
    nothing
end

function backend_change_cols_bounds!(backend::ClpJLBackend, indices::Vector{Int},
                                     lowers::Vector{Float64}, uppers::Vector{Float64})
    model = backend.model::ClpJL.ClpModel
    @inbounds for k in eachindex(indices)
        ClpJL.set_column_bounds!(model, indices[k], lowers[k], uppers[k])
    end
    nothing
end

function update_variable_block_bounds!(backend::ClpJLBackend,
                                      block::_BackendVariableBlock, offset::Int)
    model = backend.model::ClpJL.ClpModel
    @inbounds for i in eachindex(block.lower_bounds)
        ClpJL.set_column_bounds!(model, offset + i, block.lower_bounds[i], block.upper_bounds[i])
    end
    nothing
end

function update_constraint_rhs!(backend::ClpJLBackend,
                               constraint::_BackendConstraint, row::Int)
    lower, upper = _backend_constraint_bounds(constraint)
    ClpJL.set_row_bounds!(backend.model::ClpJL.ClpModel, row, lower, upper)
    nothing
end

function update_constraint_block_rhs!(backend::ClpJLBackend,
                                     block::_BackendConstraintBlock, offset::Int)
    model = backend.model::ClpJL.ClpModel
    @inbounds for (i, constraint) in enumerate(block.constraints)
        lower, upper = _backend_constraint_bounds(constraint)
        ClpJL.set_row_bounds!(model, offset + i, lower, upper)
    end
    nothing
end

function backend_update_objective_coefficients!(backend::ClpJLBackend,
                                                coefficients::_BackendObjectiveCoefficients)
    model = backend.model::ClpJL.ClpModel
    cost_factor = model.direction == ClpJL.Maximize ? -1.0 : 1.0
    @inbounds for (column, coefficient) in zip(coefficients.variable_indices,
                                              coefficients.variable_coefficients)
        scaled_coeff = cost_factor * coefficient
        model.objective[column] = scaled_coeff
        model.workspace.objective_coefficients[column] = scaled_coeff
    end
    backend.objective_changed = true
    nothing
end

function _solve_clpjl!(backend::ClpJLBackend)
    model = backend.model === nothing ?
        throw(ArgumentError("ClpJLBackend is not initialized")) : backend.model

    status = ClpJL.solve!(model; algorithm=ClpJL.DualSimplex)
    if status != ClpJL.Optimal
        status = ClpJL.solve!(model; algorithm=ClpJL.PrimalSimplex)
    end
    if status != ClpJL.Optimal
        backend.status = :infeasible
        return false, NaN
    end

    backend.objective_changed = false
    backend.objective_value = ClpJL.objective_value(model)
    backend.status = :optimal

    sol = model.workspace.primal_solution
    dual = model.workspace.dual_solution
    reduced_costs = model.workspace.reduced_costs

    num_cols = model.num_cols
    num_rows = model.num_rows

    @inbounds for i in 1:num_cols
        backend.col_value[i] = sol[i]
        backend.col_dual[i] = reduced_costs[i]
    end
    @inbounds for r in 1:num_rows
        backend.row_dual[r] = dual[r]
        backend.row_value[r] = sol[num_cols + r]
    end

    return true, backend.objective_value
end

function backend_solve!(backend::ClpJLBackend,
                        variable_blocks::Vector{_BackendVariableBlock},
                        variable_offsets::Vector{Int},
                        constraint_blocks::Vector{_BackendConstraintBlock},
                        constraint_offsets::Vector{Int})
    feasible, objective = _solve_clpjl!(backend)
    feasible || return false, objective

    for block in variable_blocks
        for i in eachindex(block.values)
            block.values[i] = backend.col_value[block.offset + i]
        end
    end
    for block in constraint_blocks
        for i in eachindex(block.values)
            block.values[i] = backend.row_dual[block.offset + i]
        end
    end

    return feasible, objective
end

end # module ReactiveModelClpJLExt
