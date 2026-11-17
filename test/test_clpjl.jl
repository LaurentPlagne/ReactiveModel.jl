using Test
using ReactiveModel
import ReactiveModel as RM
using ClpJL

@testset "ClpJLBackend: Pure-Julia Clp Integration" begin
    # 1. Economic dispatch model (Minimization)
    @model ClpJLDispatch(data, sel) begin
        G = 1:data.nG
        @variable p[g in G] in [0.0, data.cap[sel.t, g]]
        @constraint DemandBalance begin
            sum(p[g] for g in G) == data.demand[sel.t]
        end
        @objective Min begin
            sum(data.cost[sel.t, g] * p[g] for g in G)
        end
    end

    study_data = (
        nG = 2,
        cap = [80.0 60.0; 80.0 60.0],
        cost = [40.0 70.0; 40.0 70.0],
        demand = [100.0, 110.0],
    )

    pdef = ClpJLDispatch(study_data)
    cdata = RM.ConstData(pdef, study_data)

    # Test auto, dense, and sparse factorization modes
    for sparse_mode in (nothing, false, true)
        backend = ClpJLBackend(; sparse_factorization=sparse_mode)
        lp = RM.LinearProblem(cdata, (t = 1,); backend=backend)

        # Initial solve t=1: demand = 100 -> p1 = 80, p2 = 20, cost = 80*40 + 20*70 = 4600
        sol1 = RM.solve(lp, (t = 1,))
        @test RM.is_optimal(sol1)
        @test RM.objective_value(sol1) ≈ 4600.0
        @test RM.primal(sol1, :p, 1) ≈ 80.0
        @test RM.primal(sol1, :p, 2) ≈ 20.0
        @test RM.dual(sol1, :DemandBalance) ≈ 70.0
        @test RM.activity(sol1, :DemandBalance) ≈ 100.0

        # Warm solve t=2: demand = 110 -> p1 = 80, p2 = 30, cost = 80*40 + 30*70 = 5300
        sol2 = RM.solve(lp, (t = 2,))
        @test RM.is_optimal(sol2)
        @test RM.objective_value(sol2) ≈ 5300.0
        @test RM.primal(sol2, :p, 1) ≈ 80.0
        @test RM.primal(sol2, :p, 2) ≈ 30.0
        @test RM.dual(sol2, :DemandBalance) ≈ 70.0

        # Zero allocation test on warm re-solve
        allocs = @allocated RM.solve(lp, (t = 2,))
        @test allocs == 0
    end

    # 2. Parity check against Clp_jll (C-ABI)
    lp_c = RM.LinearProblem(cdata, (t = 1,); backend=ClpBackend())
    lp_jl = RM.LinearProblem(cdata, (t = 1,); backend=ClpJLBackend())

    for t_step in 1:2
        sol_c = RM.solve(lp_c, (t = t_step,))
        sol_jl = RM.solve(lp_jl, (t = t_step,))
        @test RM.objective_value(sol_c) ≈ RM.objective_value(sol_jl) atol=1e-10
        @test RM.primal(sol_c, :p, 1) ≈ RM.primal(sol_jl, :p, 1) atol=1e-10
        @test RM.primal(sol_c, :p, 2) ≈ RM.primal(sol_jl, :p, 2) atol=1e-10
        @test RM.dual(sol_c, :DemandBalance) ≈ RM.dual(sol_jl, :DemandBalance) atol=1e-10
    end

    # 3. Maximization problem with inequality constraints
    @model ClpJLMax(data, sel) begin
        @variable x1 in [0.0, 10.0]
        @variable x2 in [0.0, 10.0]
        @constraint Limit begin
            x1 + 2.0 * x2 <= sel.capacity
        end
        @objective Max begin
            3.0 * x1 + 5.0 * x2
        end
    end

    max_def = ClpJLMax(nothing)
    max_cdata = RM.ConstData(max_def, nothing)
    max_lp = RM.LinearProblem(max_cdata, (capacity = 14.0,); backend=ClpJLBackend())

    # Max 3 x1 + 5 x2 s.t. x1 + 2 x2 <= 14, 0 <= x1 <= 10, 0 <= x2 <= 10
    # Optimal: x1 = 10, x2 = 2 -> 3*10 + 5*2 = 40
    max_sol1 = RM.solve(max_lp, (capacity = 14.0,))
    @test RM.is_optimal(max_sol1)
    @test RM.objective_value(max_sol1) ≈ 40.0
    @test RM.primal(max_sol1, :x1) ≈ 10.0
    @test RM.primal(max_sol1, :x2) ≈ 2.0

    # Warm update: capacity = 20 -> x1 = 10, x2 = 5 -> 3*10 + 5*5 = 55
    max_sol2 = RM.solve(max_lp, (capacity = 20.0,))
    @test RM.is_optimal(max_sol2)
    @test RM.objective_value(max_sol2) ≈ 55.0
    @test RM.primal(max_sol2, :x1) ≈ 10.0
    @test RM.primal(max_sol2, :x2) ≈ 5.0

    # 4. Infeasible problem detection
    infeas_lp = RM.LinearProblem(max_cdata, (capacity = -5.0,); backend=ClpJLBackend())
    @test_throws RM.SolveError RM.solve(infeas_lp, (capacity = -5.0,))

    # 5. Free variable handling with Dual Simplex (Basis Crash verification)
    @model ClpJLFreeVar(data, sel) begin
        @variable x in [0.0, 10.0]
        @variable v in [-Inf, Inf]

        @constraint Cut begin
            v - 2.0 * x <= sel.rhs
        end update_rhs_when (o.rhs != n.rhs)

        @objective Max begin
            v + 3.0 * x
        end
    end

    free_def = ClpJLFreeVar(nothing)
    free_cdata = RM.ConstData(free_def, nothing)
    free_lp = RM.LinearProblem(free_cdata, (rhs = 5.0,); backend=ClpJLBackend())
    free_sol1 = RM.solve(free_lp, (rhs = 5.0,))
    @test RM.is_optimal(free_sol1)
    @test RM.objective_value(free_sol1) ≈ 55.0
    @test RM.primal(free_sol1, :x) ≈ 10.0
    @test RM.primal(free_sol1, :v) ≈ 25.0

    # Warm re-solve with Dual Simplex
    free_sol2 = RM.solve(free_lp, (rhs = 12.0,))
    @test RM.is_optimal(free_sol2)
    @test RM.objective_value(free_sol2) ≈ 62.0
    @test RM.primal(free_sol2, :x) ≈ 10.0
    @test RM.primal(free_sol2, :v) ≈ 32.0
end

@testset "ClpJL rejects unsupported matrix updates and export" begin
    @model ClpJLChangingMatrix(data, point) begin
        @variable x in [0.0, 10.0]
        @constraint floor point.a * x >= 2.0
        @objective Min x
    end
    lp = LinearProblem(ConstData(ClpJLChangingMatrix(nothing)), (a=1.0,); backend=ClpJLBackend())
    sol = solve(lp, (a=1.0,))
    @test objective_value(sol) ≈ 2.0
    @test_throws ArgumentError solve(lp, (a=2.0,))
    @test !lp.built
    @test_throws ErrorException objective_value(sol)
    rebuild!(lp)
    @test objective_value(lp, (a=2.0,)) ≈ 1.0
    mktempdir() do dir
        filename = joinpath(dir, "unsupported.lp")
        @test_throws ArgumentError backend_write_to_file(lp.backend, filename)
        @test !ispath(filename)
    end
end
