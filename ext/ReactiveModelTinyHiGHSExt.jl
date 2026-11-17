module ReactiveModelTinyHiGHSExt

using ReactiveModel
import ReactiveModel: TinyHiGHSBackend,
    _BackendVariableBlock, _BackendConstraintBlock, _BackendConstraint,
    _BackendObjective, _BackendObjectiveCoefficients, _backend_constraint_bounds, _backend_data,
    backend_has_duals, backend_set_verbose_flag!, backend_set_objective_sense!,
    backend_write_to_file, initialize!, backend_change_col_bounds!,
    backend_change_cols_bounds!, update_variable_block_bounds!,
    update_constraint_rhs!, update_constraint_variable_coefficients!,
    backend_update_objective_coefficients!, backend_solve!
import TinyHiGHS as TH

ReactiveModel.TinyHiGHSBackend() = TinyHiGHSBackend{TH.SimplexEngine,TH.WarmSolve}(
    nothing, TH.WarmSolve(), Float64[], Float64[], Float64[], Float64[], NaN, :not_solved,
)
ReactiveModel._similar_backend(::TinyHiGHSBackend) = TinyHiGHSBackend()

backend_has_duals(::TinyHiGHSBackend) = true

backend_set_verbose_flag!(::TinyHiGHSBackend; verbose::Bool) = nothing

function backend_set_objective_sense!(backend::TinyHiGHSBackend, sense::Symbol)
    highs_sense = sense === :max ? TH.kMaximize : TH.kMinimize
    TH.change_objective_sense!(backend.engine::TH.SimplexEngine, highs_sense)
    nothing
end

backend_write_to_file(backend::TinyHiGHSBackend, filename) =
    TH.write_lp(filename, (backend.engine::TH.SimplexEngine).lp)

function initialize!(backend::TinyHiGHSBackend,
                     variable_blocks::Vector{_BackendVariableBlock},
                     constraint_blocks::Vector{_BackendConstraintBlock},
                     objective::_BackendObjective)
    all(block -> block.domain === :continuous, variable_blocks) ||
        throw(ArgumentError("TinyHiGHSBackend supports continuous variables only"))
    data = _backend_data(variable_blocks, constraint_blocks, objective;
                         drop_tolerance=TH.SimplexOptions().small_matrix_value)
    ncol, nrow = data.ncols, data.nrows
    matrix_csc = data.matrix
    matrix = TH.SparseMatrix(
        ncol, nrow, matrix_csc.colptr, matrix_csc.rowval, matrix_csc.nzval,
    )
    sense = objective.op === :min ? TH.kMinimize : TH.kMaximize
    simplex_lp = TH.SimplexLp(
        ncol, nrow, matrix, data.cost, data.lower, data.upper, data.row_lower, data.row_upper;
        offset=0.0, sense=sense,
    )
    backend.engine = TH.SimplexEngine(simplex_lp, TH.SimplexOptions())
    backend.warm = TH.WarmSolve()
    backend.col_value = zeros(ncol)
    backend.col_dual = zeros(ncol)
    backend.row_value = zeros(nrow)
    backend.row_dual = zeros(nrow)
    backend.objective_value = NaN
    backend.status = :not_solved
    nothing
end

function backend_change_col_bounds!(backend::TinyHiGHSBackend, index::Int,
                                    lower::Float64, upper::Float64)
    engine = backend.engine::TH.SimplexEngine
    TH.change_col_bounds!(engine, index, lower, upper)
    nothing
end

function backend_change_cols_bounds!(backend::TinyHiGHSBackend, indices::Vector{Int},
                                     lowers::Vector{Float64}, uppers::Vector{Float64})
    TH.change_cols_bounds!(backend.engine::TH.SimplexEngine, indices, lowers, uppers)
    nothing
end

function update_variable_block_bounds!(backend::TinyHiGHSBackend,
                                      block::_BackendVariableBlock, offset::Int)
    TH.change_cols_bounds!(backend.engine::TH.SimplexEngine,
                           (offset + 1):(offset + length(block.lower_bounds)),
                           block.lower_bounds, block.upper_bounds)
end

function update_constraint_rhs!(backend::TinyHiGHSBackend,
                               constraint::_BackendConstraint, row::Int)
    lower, upper = _backend_constraint_bounds(constraint)
    TH.change_row_bounds!(backend.engine::TH.SimplexEngine, row, lower, upper)
    nothing
end

function update_constraint_variable_coefficients!(backend::TinyHiGHSBackend,
                                                  row::Int,
                                                  constraint::_BackendConstraint)
    engine = backend.engine::TH.SimplexEngine
    for (column, coefficient) in zip(constraint.variable_indices,
                                     constraint.variable_coefficients)
        TH.change_coeff!(engine, row, column, coefficient)
    end
    nothing
end

function backend_update_objective_coefficients!(backend::TinyHiGHSBackend,
                                                coefficients::_BackendObjectiveCoefficients)
    TH.change_cols_cost!(
        backend.engine::TH.SimplexEngine,
        coefficients.variable_indices, coefficients.variable_coefficients,
    )
    nothing
end

function _solve_tiny!(backend::TinyHiGHSBackend)
    engine = backend.engine === nothing ?
        throw(ArgumentError("TinyHiGHSBackend is not initialized")) : backend.engine
    # A warm basis that ends a solve without a verdict restarts from the last
    # optimal basis, then the logical one (see `TH.warm_solve!`).
    TH.warm_solve!(backend.warm, engine)
    TH.restore_scale!(engine)
    if engine.model_status != TH.kOptimal
        backend.status = :infeasible
        return false, NaN
    end
    TH.compute_primal_objective_value!(engine)
    sense = Int(engine.lp.sense)
    @inbounds for i in eachindex(backend.col_value)
        backend.col_value[i] = engine.info.workValue[i]
        backend.col_dual[i] = -sense * engine.info.workDual[i]
    end
    @inbounds for row in 1:engine.lp.num_row
        variable = engine.basis.basicIndex[row]
        if variable <= engine.lp.num_col
            backend.col_value[variable] = engine.info.baseValue[row]
            backend.row_value[row] = -engine.info.workValue[engine.lp.num_col + row]
        else
            backend.row_value[row] = -engine.info.baseValue[row]
        end
        backend.row_dual[row] = -sense * engine.info.workDual[engine.lp.num_col + row]
    end
    backend.objective_value = engine.info.primal_objective_value
    backend.status = :optimal
    true, backend.objective_value
end

function backend_solve!(backend::TinyHiGHSBackend,
                        variable_blocks::Vector{_BackendVariableBlock},
                        variable_offsets::Vector{Int},
                        constraint_blocks::Vector{_BackendConstraintBlock},
                        constraint_offsets::Vector{Int})
    feasible, objective = _solve_tiny!(backend)
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
    feasible, objective
end

end # module
