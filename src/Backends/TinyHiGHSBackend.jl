"""
    TinyHiGHSBackend()

Optional persistent pure-Julia TinyHiGHS backend for continuous linear programs.
Install TinyHiGHS separately and load it with `using TinyHiGHS` before constructing
this backend. HiGHSBackend and ClpBackend are available without optional solvers.
"""
mutable struct TinyHiGHSBackend{E,W} <: AbstractBackend
    engine::Union{Nothing,E}
    warm::W
    col_value::Vector{Float64}
    col_dual::Vector{Float64}
    row_value::Vector{Float64}
    row_dual::Vector{Float64}
    objective_value::Float64
    status::Symbol
end

function TinyHiGHSBackend(args...; kwargs...)
    throw(ArgumentError("TinyHiGHSBackend requires TinyHiGHS.jl; install it separately and run `using TinyHiGHS` first"))
end
