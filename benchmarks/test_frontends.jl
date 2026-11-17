"""
Frontend regression checks on independently specified synthetic cascades.

    julia --startup-file=no --project=benchmarks benchmarks/test_frontends.jl [preset ...]

Check fixed public preset dimensions, numerical parity and warmed allocation
ceilings. Times are reported rather than asserted. Allocation bounds apply to
these implementations and dependency versions, not to all client applications.
"""

using Test

include(joinpath(@__DIR__, "valley_comparison.jl"))

# Public synthetic experiment sizes: columns, rows and sparse entries.
const SYNTHETIC_SIZES = Dict(
    "ValleySmall" => (ncols=50, nrows=44),
    "ValleyMedium" => (ncols=172, nrows=102),
    "ValleyLarge" => (ncols=390, nrows=236),
)

# Regression ceilings with margin, not claims about other implementations.
const ALLOCATION_CEILINGS = Dict(
    :rm_highs => 1_024,
    :rm_clp => 1_024,
    :moi => 64_000,
    :jump => 128_000,
    :jump_user_pred => 256_000,
    :jump_user => 512_000,
)

"""
    check_sizes(data)

The generated LP has the documented public experiment dimensions.
"""
function check_sizes(data::ValleyData)
    described = SV.describe(data)
    expected = SYNTHETIC_SIZES[described.name]
    @test described.ncols == expected.ncols
    @test described.nrows == expected.nrows
    described
end

"""
    check_parity(models, points)

Every front-end returns ReactiveModel's objective and its first water-balance
dual, at every point.

The dual is what makes this check bite. A front-end whose updater misses a
moving right-hand side still returns a plausible objective; the dual of the row
that did not move is what exposes it. JuMP, MOI and POI report MOI's opposite
sign convention, hence the negation.
"""
function check_parity(models, points; atol=1e-7, rtol=1e-9)
    row = models.l.water_balance[1, 1]
    worst_objective, worst_dual = 0.0, 0.0
    for point in points
        solution = RM.solve(models.rm_highs, point)
        reference = RM.objective_value(solution)
        reference_dual = RM.dual(solution, :water_balance, (1, 1))

        objectives = (
            jump=jump_solve!(models.jump, point),
            jump_user=jump_user_solve!(models.jump_user, point),
            jump_user_pred=jump_user_solve!(models.jump_user_pred, point;
                predicated=true),
            moi=moi_solve!(models.moi, point),
            poi=poi_solve!(models.poi, point),
            poi_hybrid=poi_solve!(models.poi_hybrid, point),
        )
        for (name, value) in pairs(objectives)
            isapprox(reference, value; atol, rtol) ||
                error("objective mismatch on $name at $point: " *
                      "$reference versus $value")
            worst_objective = max(worst_objective, abs(reference - value))
        end

        duals = (
            jump=JuMP.dual(models.jump.rowrefs[row]),
            jump_user=JuMP.dual(models.jump_user.water_balance[1, 1]),
            jump_user_pred=JuMP.dual(models.jump_user_pred.water_balance[1, 1]),
            moi=MOI.get(models.moi.optimizer, MOI.ConstraintDual(),
                models.moi.rowrefs[row]),
            poi=JuMP.dual(models.poi.rowrefs[row]),
            poi_hybrid=JuMP.dual(models.poi_hybrid.rowrefs[row]),
        )
        for (name, value) in pairs(duals)
            isapprox(-reference_dual, value; atol=1e-6, rtol=1e-6) ||
                error("dual mismatch on $name at $point: " *
                      "$(-reference_dual) versus $value")
            worst_dual = max(worst_dual, abs(-reference_dual - value))
        end
    end
    (; points=length(points), worst_objective, worst_dual)
end

"""
    check_tiny_parity(models, points)

TinyHiGHS is a different simplex, so it is held to a looser tolerance. At a
degenerate optimum it may also select a different valid dual, so only the
objective is asserted.
"""
function check_tiny_parity(models, points; tiny_atol=1e-5)
    worst = 0.0
    for point in points
        reference = RM.objective_value(models.rm_highs, point)
        worst = max(worst,
            abs(reference - RM.objective_value(models.rm_tiny, point)))
    end
    @test worst <= tiny_atol * max(1.0, worst)
    worst
end

"""
    check_clp_parity(models, points)

`ClpBackend` is a second implementation of the same direct-backend contract
`HiGHSBackend` implements, not a modeling layer, so it is checked the same way
`rm_tiny` is: against ReactiveModel's own HiGHS solve, not against JuMP/MOI/POI.
Held to a tighter tolerance than TinyHiGHS since Clp and HiGHS are both
full-precision revised-simplex implementations; the dual is reported rather
than asserted, since a degenerate optimum can have more than one valid dual and
two different simplex codes are not guaranteed to land on the same one.
"""
function check_clp_parity(models, points; clp_atol=1e-6)
    worst_objective = 0.0
    worst_dual = 0.0
    row_key = (1, 1)
    for point in points
        reference = RM.solve(models.rm_highs, point)
        clp = RM.solve(models.rm_clp, point)
        worst_objective = max(worst_objective,
            abs(RM.objective_value(reference) - RM.objective_value(clp)))
        worst_dual = max(worst_dual,
            abs(RM.dual(reference, :water_balance, row_key) -
                RM.dual(clp, :water_balance, row_key)))
    end
    @test worst_objective <= clp_atol * max(1.0, worst_objective)
    (; worst_objective, worst_dual)
end

"""
    measure(models, path, points; passes)

Per-LP bytes and microseconds over a warmed steady-state loop. The first pass is
outside the measurement, so no path pays for its cold basis or its compilation.
"""
function measure(models, path::Symbol, points; passes::Int=3)
    run_path!(models, path, points; read_dual=false)
    GC.gc()
    sample = @timed (for _ in 1:passes
        run_path!(models, path, points; read_dual=false)
    end)
    n = length(points) * passes
    (; path, bytes=sample.bytes / n, us=1e6 * sample.time / n)
end

function test_preset(spec::SV.ValleySpec; nscenarios::Int=6,
        stocks_per_scenario::Int=5)
    data = SV.valley_data(spec; nregions=1)
    described = check_sizes(data)
    points = SV.valley_points(data; nscenarios, order=:refinement,
        stocks_per_scenario)

    models = build_models(data, points)

    parity = check_parity(models, points)
    @test parity.points == length(points)
    tiny_error = BENCH_TINY ? check_tiny_parity(models, points) : nothing
    clp_parity = check_clp_parity(models, points)

    measurements = [measure(models, path, points)
                    for path in keys(ALLOCATION_CEILINGS)]
    by_path = Dict(m.path => m for m in measurements)

    for (path, ceiling) in ALLOCATION_CEILINGS
        @test by_path[path].bytes <= ceiling
    end

    # The ordering itself, which is what the parallel-scaling conclusion rests
    # on: ReactiveModel below direct MOI, direct MOI below matrix-driven JuMP,
    # and the naive user-written update above the predicated one.
    @test by_path[:rm_highs].bytes < by_path[:moi].bytes
    @test by_path[:moi].bytes < by_path[:jump].bytes
    @test by_path[:jump_user_pred].bytes < by_path[:jump_user].bytes

    println("  ", described.name, ": ", described.ncols, " cols / ",
        described.nrows, " rows / ", described.nnz, " nnz, ",
        parity.points, " points")
    println("    parity: objective <= ", round(parity.worst_objective; sigdigits=2),
        ", dual <= ", round(parity.worst_dual; sigdigits=2),
        ", TinyHiGHS objective <= ", (tiny_error === nothing ? "disabled" : round(tiny_error; sigdigits=2)),
        ", Clp objective <= ", round(clp_parity.worst_objective; sigdigits=2),
        ", Clp dual <= ", round(clp_parity.worst_dual; sigdigits=2))
    for m in sort(measurements; by=m -> m.bytes)
        println("    ", rpad(String(m.path), 16),
            lpad(round(Int, m.bytes), 8), " B/LP  ",
            lpad(round(m.us; digits=1), 8), " us/LP")
    end
    (; described, parity, measurements)
end

"""
    test_jump_expressiveness()

Two affine formulations of the synthetic transition LP and JuMP macro syntax
will not spell, pinned as tests rather than left as prose in a report.

Neither is a defect: both follow from reasonable JuMP design decisions. They are
recorded because a reader comparing the two front-ends should know that the
workarounds in `jump_user_model.jl` are necessary, and because a future JuMP
release that removes the need for them should show up as a failing test here,
not as a footnote nobody rechecks.
"""
function test_jump_expressiveness()
    function_type(reference) =
        typeof(JuMP.index(reference)).parameters[1]
    is_row(reference) =
        function_type(reference) === MOI.ScalarAffineFunction{Float64}
    is_bound(reference) = function_type(reference) === MOI.VariableIndex

    model = JuMP.Model()
    JuMP.@variable(model, 0 <= x <= 10)

    # Step 1. An infinite right-hand side is refused by the operator form. The
    # macro moves the bound into the affine expression's constant, and an
    # infinite constant makes it `NaN`. The synthetic LP keeps rows neutralized by an
    # infinite right-hand side so the row count and the pivoting stay stable
    # across points, so this is not an exotic need.
    @test_throws Exception JuMP.@constraint(model, x <= Inf)
    @test_throws Exception JuMP.@constraint(model, 1.0 * x <= Inf)

    # Step 2. That pushes the user to the set form — which, on a *bare
    # variable*, JuMP turns into a variable bound rather than a row. Here it
    # collides with the upper bound `x` already carries; on a variable with no
    # bound it would silently become a bound with no dual and no place in the
    # matrix, which is not what a neutralized row is.
    @test_throws MOI.UpperBoundAlreadySet JuMP.@constraint(
        model, x in MOI.LessThan(Inf))

    unbounded = JuMP.Model()
    JuMP.@variable(unbounded, y)
    @test is_bound(JuMP.@constraint(unbounded, y in MOI.LessThan(4.0)))

    # Step 3. Forcing the variable into an affine expression is the way out, and
    # is all `as_row` does.
    neutralized = JuMP.@constraint(model, as_row(x) in MOI.LessThan(Inf))
    @test is_row(neutralized)
    @test JuMP.normalized_rhs(neutralized) == Inf

    # For the record, and against the wording of earlier reports: with a
    # *finite* right-hand side the operator form does produce a row on a bare
    # variable in this JuMP version, bound or no bound. The trap above is
    # specific to the set form, which the infinite right-hand side forces.
    @test is_row(JuMP.@constraint(model, x <= 4.0))
    @test is_row(JuMP.@constraint(model, x >= 1.0))
    return nothing
end

const TEST_PRESETS = Dict(
    "small" => SV.VALLEY_SMALL,
    "medium" => SV.VALLEY_MEDIUM,
    "large" => SV.VALLEY_LARGE,
)

function main(names=isempty(ARGS) ? ["small", "medium", "large"] : ARGS)
    println((julia=VERSION, cpu=Sys.cpu_info()[1].model,
        threads=Threads.nthreads()))
@testset "independent synthetic frontends" begin
        @testset "what JuMP will not express" begin
            test_jump_expressiveness()
        end
        for name in names
            @testset "$name" begin
                test_preset(TEST_PRESETS[lowercase(name)])
            end
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
