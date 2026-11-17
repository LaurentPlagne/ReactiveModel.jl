# ---------------------------------------------------------------------------
# XpressBackend: FICO Xpress backend for MILP and LP
# ---------------------------------------------------------------------------

"""
    XpressBackend(; kwargs...)

Persistent FICO Xpress backend for Mixed-Integer Linear Programs (MILP) and LPs,
wrapping `Xpress.jl` native C API (`XPRSloadmip`, `XPRSmipoptimize`, etc.).

Requires `using Xpress` to load the `ReactiveModelXpressExt` extension.
This optional integration has not been qualified against a licensed solver.
"""
mutable struct XpressBackend <: AbstractBackend
    prob::Any
    col_value::Vector{Float64}
    col_dual::Vector{Float64}
    row_value::Vector{Float64}
    row_dual::Vector{Float64}
    a_start::Vector{Cint}
    a_index::Vector{Cint}
    a_value::Vector{Cdouble}
    col_cost::Vector{Cdouble}
    col_lower::Vector{Cdouble}
    col_upper::Vector{Cdouble}
    row_lower::Vector{Cdouble}
    row_upper::Vector{Cdouble}
    rowtype::Vector{UInt8}
    rhs::Vector{Cdouble}
    integrality::Vector{UInt8}
    int_indices::Vector{Cint}
    is_mip::Bool
    objective_value::Float64
    status::Symbol
    options::Dict{Symbol, Any}

    function XpressBackend(prob, col_value, col_dual, row_value, row_dual,
                           a_start, a_index, a_value, col_cost, col_lower, col_upper,
                           row_lower, row_upper, rowtype, rhs, integrality, int_indices,
                           is_mip, objective_value, status, options)
        new(prob, col_value, col_dual, row_value, row_dual,
            a_start, a_index, a_value, col_cost, col_lower, col_upper,
            row_lower, row_upper, rowtype, rhs, integrality, int_indices,
            is_mip, objective_value, status, options)
    end
end

function XpressBackend(args...; kwargs...)
    error("`XpressBackend` requires `Xpress.jl` to be loaded.\nPlease run `using Xpress` or `import Xpress` first.")
end

# Status codes from Xpress.jl's generated C API. A nonzero stop reason is
# distinct from a proof of infeasibility. MIP_UNBOUNDED describes an unbounded
# continuous relaxation, so it does not establish that the MILP is unbounded.
function _xpress_status(is_mip::Bool, status::Integer, stop::Integer)
    stop == 1 && return :time_limit
    stop in (2, 9) && return :interrupted
    stop == 3 && return :node_limit
    stop == 4 && return :iteration_limit
    stop == 6 && return :solution_limit
    stop == 8 && return :memory_limit
    stop == 13 && return :numerical_error
    stop == 14 && return :work_limit
    stop in (0, 5, 10) || return :other_error
    if is_mip
        status == 6 && return :optimal
        status == 5 && return :infeasible
        status == 7 && return :infeasible_or_unbounded
        status == 4 && return :feasible
    else
        status == 1 && return :optimal
        status == 2 && return :infeasible
        status == 5 && return :unbounded
        status in (3, 6) && return :objective_limit
        status == 7 && return :numerical_error
    end
    return :other_error
end
