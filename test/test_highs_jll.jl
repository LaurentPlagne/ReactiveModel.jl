using Test
using ReactiveModel
import ReactiveModel as RM
using HiGHS_jll

@testset "HiGHS_jll: Direct C-ABI and Backend Parity" begin
    # 1. Direct C-ABI call verification via HiGHS_jll
    highs = ccall((:Highs_create, HiGHS_jll.libhighs), Ptr{Cvoid}, ())
    @test highs != C_NULL
    status = ccall((:Highs_setStringOptionValue, HiGHS_jll.libhighs), Cint, (Ptr{Cvoid}, Cstring, Cstring), highs, "output_flag", "false")
    @test status == 0
    ccall((:Highs_destroy, HiGHS_jll.libhighs), Cvoid, (Ptr{Cvoid},), highs)

    # 2. Multi-backend solve parity: direct HiGHS_jll vs Clp_jll.
    # The optional ClpJL integration has its own test file.
    @model TestDispatchJLL(data, sel) begin
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

    pdef = TestDispatchJLL(study_data)
    cdata = RM.ConstData(pdef, study_data)

    lp_highs = RM.LinearProblem(cdata, (t = 1,); backend=HiGHSBackend())
    lp_clp = RM.LinearProblem(cdata, (t = 1,); backend=ClpBackend())

    for t_step in 1:2
        sol_h = RM.solve(lp_highs, (t = t_step,))
        sol_c = RM.solve(lp_clp, (t = t_step,))

        @test RM.is_optimal(sol_h)
        @test RM.is_optimal(sol_c)

        @test RM.objective_value(sol_h) ≈ RM.objective_value(sol_c) atol=1e-10

        @test RM.primal(sol_h, :p, 1) ≈ RM.primal(sol_c, :p, 1) atol=1e-10
        @test RM.primal(sol_h, :p, 2) ≈ RM.primal(sol_c, :p, 2) atol=1e-10

        @test RM.dual(sol_h, :DemandBalance) ≈ RM.dual(sol_c, :DemandBalance) atol=1e-10
    end

    # 3. Warm re-solve performance with HiGHSBackend
    sol_warm = RM.solve(lp_highs, (t = 2,))
    allocs = @allocated RM.solve(lp_highs, (t = 2,))
    @test allocs == 0
end

@testset "HiGHS_jll: options, integrality, updates and resource lifecycle" begin
    backend = HiGHSBackend(output_flag=false, simplex_strategy=1,
                           time_limit=10.0, presolve="off")
    @test RM._HiGHSCAPI.Highs_getSizeofHighsInt(backend) == sizeof(Cint)
    @test_throws ErrorException HiGHSBackend(not_a_highs_option=1)
    compatible_options = HiGHSBackend(time_limit=10, output_flag="false", presolve=:off)
    time_limit = Ref{Cdouble}(NaN)
    @test ccall((:Highs_getDoubleOptionValue, HiGHS_jll.libhighs), Cint,
        (Ptr{Cvoid}, Cstring, Ptr{Cdouble}), compatible_options, "time_limit", time_limit) == 0
    @test time_limit[] == 10.0
    finalize(compatible_options)

    # Continuous relaxation gives 1.5; native integer and binary models give 2.
    for domain in (:integer, :binary)
        definition = ProblemDefinition(
            variables=(x=VariableBlock(:x, 1:3; lower=0.0, upper=1.0, domain),),
            constraints=(floor=ConstraintBlock(:floor, (nothing,), :ge,
                (eq, d, p, k) -> begin
                    for i in 1:3
                        set_coeff!(eq, :x, i, 1.0)
                    end
                end, (d, p, k) -> p.floor),),
            objective=Objective(:min, (eq, d, p) -> begin
                for i in 1:3
                    set_coeff!(eq, :x, i, 1.0)
                end
            end),
        )
        lp = LinearProblem(definition, (floor=1.5,); backend=HiGHSBackend())
        @test objective_value(lp, (floor=1.5,)) ≈ 2.0
        @test objective_value(lp, (floor=2.5,)) ≈ 3.0
        @test !backend_has_duals(lp.backend)
        @test all(isnan, lp.backend.row_dual)
    end

    # More than three columns and rows exercise the C range-based update paths;
    # objective and matrix coefficients also change between the two solves.
    definition = ProblemDefinition(
        variables=(x=VariableBlock(:x, 1:4; lower=0.0,
            upper=(d, p, i) -> p.upper),),
        constraints=(cap=ConstraintBlock(:cap, 1:4, :le,
            (eq, d, p, i) -> set_coeff!(eq, :x, i, p.coefficient),
            (d, p, i) -> p.rhs),),
        objective=Objective(:max, (eq, d, p) -> begin
            for i in 1:4
                set_coeff!(eq, :x, i, p.cost)
            end
        end),
    )
    p1 = (upper=5.0, coefficient=1.0, rhs=4.0, cost=1.0)
    p2 = (upper=3.0, coefficient=2.0, rhs=5.0, cost=2.0)
    lp = LinearProblem(definition, p1; backend)
    @test objective_value(lp, p1) ≈ 16.0
    @test objective_value(lp, p2) ≈ 20.0
    @test all(i -> primal(solve(lp, p2), :x, i) ≈ 2.5, 1:4)
    backend_change_cols_bounds!(backend, [1, 3], [0.0, 0.0], [2.0, 2.0])
    @test first(backend_solve!(backend, lp.backend_variable_blocks,
        lp.variable_block_offsets, lp.backend_constraint_blocks, lp.constraint_block_offsets))
    @test backend.objective_value ≈ 18.0
    backend_set_objective_sense!(backend, :min)
    @test first(backend_solve!(backend, lp.backend_variable_blocks,
        lp.variable_block_offsets, lp.backend_constraint_blocks, lp.constraint_block_offsets))
    @test backend.objective_value ≈ 0.0

    mktempdir() do dir
        filename = joinpath(dir, "direct-highs.lp")
        # HiGHS can return a warning when generating default row/column names.
        @test backend_write_to_file(backend, filename) >= 0
        @test isfile(filename)
    end
    finalize(backend)
    @test backend.inner == C_NULL
    RM._highs_finalize!(backend) # The cleanup operation is idempotent.
    @test backend.inner == C_NULL
end
