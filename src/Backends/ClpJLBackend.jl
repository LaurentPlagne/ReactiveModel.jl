# ---------------------------------------------------------------------------
# ClpJLBackend: Pure-Julia Simplex backend for continuous linear programs
# ---------------------------------------------------------------------------

"""
    ClpJLBackend(; sparse_factorization=nothing, options=...)

Persistent pure-Julia Clp backend for continuous linear programs, wrapping `ClpJL.jl`.
It eliminates C-FFI call boundaries (`ccall`) and executes directly in native Julia threads.
When `sparse_factorization` is `nothing` (default), sparse Markowitz factorization with Forrest-Tomlin
updates is automatically selected for models with more than 40 constraints (`nrow > 40`), and dense
factorization otherwise.

Matrix coefficient updates and model export are not supported and raise an
`ArgumentError`. Use HiGHSBackend or ClpBackend for these operations, or call
`rebuild!` before solving after a matrix change.

Requires `using ClpJL` to load the `ReactiveModelClpJLExt` extension.
"""
mutable struct ClpJLBackend <: AbstractBackend
    model::Any
    sparse_factorization::Union{Nothing, Bool}
    options::Any
    objective_changed::Bool
    col_value::Vector{Float64}
    col_dual::Vector{Float64}
    row_value::Vector{Float64}
    row_dual::Vector{Float64}
    objective_value::Float64
    status::Symbol

    function ClpJLBackend(model, sparse_factorization, options, objective_changed,
                          col_value, col_dual, row_value, row_dual, objective_value, status)
        new(model, sparse_factorization, options, objective_changed,
            col_value, col_dual, row_value, row_dual, objective_value, status)
    end
end

function ClpJLBackend(args...; kwargs...)
    error("`ClpJLBackend` requires `ClpJL.jl` to be loaded.\nPlease run `using ClpJL` or `import ClpJL` first.")
end

function update_constraint_variable_coefficients!(::ClpJLBackend, ::Int, ::_BackendConstraint)
    throw(ArgumentError("ClpJLBackend does not support matrix coefficient updates; use HiGHSBackend or ClpBackend, or call rebuild! before solving"))
end

function backend_write_to_file(::ClpJLBackend, filename)
    throw(ArgumentError("ClpJLBackend does not support model export; use HiGHSBackend or ClpBackend"))
end
