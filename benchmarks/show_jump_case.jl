"""
Run the measured JuMP formulations with HiGHS and Clp.

    julia --startup-file=no --project=benchmarks benchmarks/show_jump_case.jl [highs|clp]

No argument checks both solvers. Flat and named/block-guarded models are built
once, then updated at six points and compared with the corresponding direct RM
JLL backend. The Clp case uses the campaign's dual/no-presolve options.
"""
module JuMPBenchmarkWalkthrough
include("parallel_comparison.jl")
const P = ParallelSolverComparison

function main(; solvers=(:highs, :clp))
    data = P.SV.valley_data(P.SV.VALLEY_SMALL; nregions=1)
    points = P.SV.valley_points(data; nscenarios=2, order=:refinement, stocks_per_scenario=3)
    layout = P.SV.layout(data)
    println("Synthetic reservoir LP: ", layout.ncols, " variables, ", layout.nrows,
            " rows, ", data.spec.npostes, " periods")
    for solver in solvers
        solver in (:highs, :clp) || throw(ArgumentError("choose highs or clp"))
        flat_tag = solver === :clp ? Val(:jump_clp_flat) : Val(:jump_flat)
        named_tag = solver === :clp ? Val(:jump_clp_named_guarded) : Val(:jump_named_guarded)
        rm_tag = solver === :clp ? Val(:rm_clp) : Val(:rm_highs)
        flat = P.make_path(flat_tag, data, layout, first(points))
        named = P.make_path(named_tag, data, layout, first(points))
        rm = P.make_path(rm_tag, data, layout, first(points))
        println(solver, ": point | flat JuMP objective | named JuMP objective")
        for (i, point) in enumerate(points)
            flat_value = P.solve_point!(flat, point)
            named_value = P.solve_point!(named, point)
            reference = P.solve_point!(rm, point)
            @assert isapprox(flat_value, reference; atol=1e-6, rtol=1e-8)
            @assert isapprox(named_value, reference; atol=1e-6, rtol=1e-8)
            println(i, " | ", flat_value, " | ", named_value)
        end
    end
    nothing
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    solvers = isempty(ARGS) ? (:highs, :clp) : (Symbol(only(ARGS)),)
    JuMPBenchmarkWalkthrough.main(; solvers)
end
