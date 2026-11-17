using Test
using ReactiveModel
import ReactiveModel: lower_bound, upper_bound, lhs_terms!, rhs_value, objective_terms!, @coeff
include(joinpath(@__DIR__, "..", "examples", "examples.jl"))
using .ReactiveModelExamples

struct BinaryDecision <: AbstractVariable end
ReactiveModel.lower_bound(::BinaryDecision, data, parameters, _) = 0.0
ReactiveModel.upper_bound(::BinaryDecision, data, parameters, _) = 1.0

struct BinaryFloor <: AbstractConstraint end
ReactiveModel.lhs_terms!(::BinaryFloor, equation, data, parameters, _) =
    (@coeff equation decision ← 1.0)
ReactiveModel.rhs_value(::BinaryFloor, data, parameters, _) = 0.5

struct BinaryCost <: AbstractObjective end
ReactiveModel.objective_terms!(::BinaryCost, expression, data, parameters) =
    (@coeff expression decision ← 1.0)

struct MissingUpperBound <: AbstractVariable end
ReactiveModel.lower_bound(::MissingUpperBound, data, parameters, _) = 0.0

struct ChangingSupportPoint
    extra::Bool
end

struct ChangingSupportVariable <: AbstractVariable end
ReactiveModel.lower_bound(::ChangingSupportVariable, data, parameters, _) = 0.0
ReactiveModel.upper_bound(::ChangingSupportVariable, data, parameters, _) = 1.0

struct ChangingSupportCap <: AbstractConstraint end
function ReactiveModel.lhs_terms!(::ChangingSupportCap, equation, data,
                                  point::ChangingSupportPoint, _)
    @coeff equation x ← 1.0
    point.extra && (@coeff equation y ← 1.0)
end
ReactiveModel.rhs_value(::ChangingSupportCap, data, parameters, _) = 1.0

struct ChangingSupportGain <: AbstractObjective end
function ReactiveModel.objective_terms!(::ChangingSupportGain, expression, data,
                                        point::ChangingSupportPoint)
    @coeff expression x ← 1.0
    point.extra && (@coeff expression y ← 2.0)
end

struct BatchRowPoint
    step::Int
end

struct IrrelevantPoint
    id::Int
end

struct BatchRow <: AbstractConstraint
    lhs_calls::Base.RefValue{Int}
end
function ReactiveModel.lhs_terms!(row::BatchRow, equation, data, point, key)
    row.lhs_calls[] += 1
    @coeff equation x[key] ← 1.0
end
ReactiveModel.rhs_value(::BatchRow, data, point::BatchRowPoint, key) =
    Float64(key + point.step)
ReactiveModel.lhs_has_changed(::BatchRow, old::BatchRowPoint, new::BatchRowPoint) = false
ReactiveModel.rhs_has_changed(::BatchRow, old::BatchRowPoint, new::BatchRowPoint) =
    old.step != new.step

struct BatchRowObjective <: AbstractObjective end
function ReactiveModel.objective_terms!(::BatchRowObjective, expression, data, point)
    for key in 1:4
        @coeff expression x[key] ← 1.0
    end
end

dispatch_study() = dispatch_const_data(
    [100.0, 110.0, 110.0, 200.0],
    [80.0 60.0; 80.0 60.0; 80.0 60.0; 80.0 60.0],
    [40.0 70.0; 40.0 70.0; 45.0 70.0; 40.0 70.0],
)

@testset "Prebuilt model solves after a numerically unchanged parameter update" begin
    definition = ProblemDefinition(
        variables=(x=VariableBlock(:x; lower=0.0, upper=10.0),),
        constraints=(cap=ConstraintBlock(
            :cap, (nothing,), :le,
            (emit, data, point, key) -> (@coeff emit x ← 1.0),
            (data, point, key) -> 7.0,
        ),),
        objective=Objective(:max, (emit, data, point) -> (@coeff emit x ← 1.0)),
    )
    first_point, next_point = IrrelevantPoint(1), IrrelevantPoint(2)
    model = LinearProblem(ConstData(definition), first_point;
        backend=HiGHSBackend(output_flag=false))
    ReactiveModel._build!(model, first_point)

    @test model.last_status === :not_solved
    @test objective_value(solve(model, next_point)) == 7.0
    @test primal(solve(model, next_point), :x) == 7.0
end

@testset "Missing mathematical methods fail explicitly" begin
    definition = ProblemDefinition(
        variables=(decision=VariableBlock(MissingUpperBound),),
        constraints=(floor=ConstraintBlock(:ge, BinaryFloor),),
        objective=Objective(:min, BinaryCost),
    )
    model = LinearProblem(
        ConstData(definition),
        nothing;
        backend=HiGHSBackend(output_flag=false),
    )
    @test_throws MethodError solve(model, nothing)
end

@testset "Typed empty blocks are structural holes" begin
    definition = ProblemDefinition(
        variables=(
            x=VariableBlock(:x; lower=0.0, upper=10.0),
            hole=VariableBlock(
                :hole,
                Int[];
                lower=(data, parameters, key) -> error("empty variable block evaluated"),
                upper=(data, parameters, key) -> error("empty variable block evaluated"),
            ),
        ),
        constraints=(
            cap=ConstraintBlock(
                :cap,
                (nothing,),
                :le,
                (emit, data, parameters, key) -> (@coeff emit x ← 1.0),
                (data, parameters, key) -> 7.0,
            ),
            hole=ConstraintBlock(
                :hole,
                Tuple{Int}[],
                :eq,
                (emit, data, parameters, key) -> error("empty constraint block evaluated"),
                (data, parameters, key) -> error("empty constraint block evaluated"),
            ),
        ),
        objective=Objective(:max, (emit, data, parameters) -> (@coeff emit x ← 1.0)),
    )

    @test isempty(definition.variables.hole.keys)
    @test eltype(definition.variables.hole.keys) === Int
    @test isempty(definition.constraints.hole.keys)
    @test eltype(definition.constraints.hole.keys) === Tuple{Int}
    @test_throws ArgumentError VariableBlock(:untyped, Any[]; lower=0.0, upper=1.0)
    @test_throws ArgumentError ConstraintBlock(
        :untyped,
        Any[],
        :eq,
        (emit, data, parameters, key) -> nothing,
        (data, parameters, key) -> 0.0,
    )

    for backend in (HiGHSBackend(output_flag=false), ClpBackend())
        solution = solve(LinearProblem(ConstData(definition), nothing; backend), nothing)
        @test objective_value(solution) == 7.0
        @test primal(solution, :x) == 7.0
        @test isempty(solution.variables.hole.values)
        @test isempty(solution.constraints.hole.duals)
        @test isempty(solution.constraints.hole.activities)
    end
end

@testset "ReactiveModel data separation" begin
    p = DispatchParameters(1)
    const_data = dispatch_study()
    lp = LinearProblem(const_data, Parameters(p); backend=HiGHSBackend(output_flag=false))

    @test const_data isa ConstData
    @test const_data.data.generator_indices == 1:2
    @test const_data.data.demand == [100.0, 110.0, 110.0, 200.0]
    @test const_data.data.capacity isa Matrix{Float64}
    @test const_data.data.cost isa Matrix{Float64}
    @test lp.const_data === const_data
    @test lp.workspace isa Workspace
    @test lp.const_data.definition.variables.production.variable isa AbstractVariable
    @test :variable ∈ fieldnames(VariableBlock)
    @test eltype(lp.const_data.definition.variables.production.keys) === Int
    @test isconcretetype(fieldtype(typeof(lp.const_data.definition), :variables))
    @test fieldnames(DispatchParameters) == (:time_step,)
    @test Parameters(p).time_step == 1
    @test isbitstype(typeof(p))
    @test isbitstype(typeof(Parameters(p)))
    @test_throws ArgumentError Parameters([1])

    compact_definition = ProblemDefinition(
        variables=(production=VariableBlock(
            ReactiveModelExamples.DispatchProduction, 1:2),),
        constraints=(balance=ConstraintBlock(
            (nothing,), :eq, ReactiveModelExamples.DispatchBalance),),
        objective=Objective(:min, ReactiveModelExamples.DispatchObjective),
    )
    @test compact_definition.variables.production.name === :production
    @test compact_definition.constraints.balance.name === :balance
    @test length(compact_definition.constraints.balance.keys) == 1
    @test only(compact_definition.constraints.balance.keys) === nothing
    @test compact_definition.objective.objective isa
        ReactiveModelExamples.DispatchObjective
    @test_throws ArgumentError ProblemDefinition(
        variables=(production=VariableBlock(
            :wrong_name, 1:2, ReactiveModelExamples.DispatchProduction),),
        constraints=(balance=ConstraintBlock(
            (nothing,), :eq, ReactiveModelExamples.DispatchBalance),),
        objective=Objective(:min, ReactiveModelExamples.DispatchObjective),
    )

    first_solution = @inferred solve(lp, p)
    @test first_solution.variables.production[1] == 80.0
    @test first_solution.variables.production[2] == 20.0
    @test objective_value(first_solution) == 4600.0
    @test solve(lp, DispatchParameters(1)) === first_solution

    changed = solve(lp, DispatchParameters(2))
    @test changed.variables.production[2] == 30.0
    @test objective_value(changed) == 5300.0
    @test solve(lp, DispatchParameters(2)) === changed

    new_costs = solve(lp, DispatchParameters(3))
    @test new_costs.variables.production[1] == 80.0
    @test objective_value(new_costs) == 5700.0

    closure_definition = ProblemDefinition(
        variables=(production=VariableBlock(:production, 1:2;
            lower=0.0, upper=(c, parameters, g) -> c.capacity[parameters.time_step, g]),),
        constraints=(balance=ConstraintBlock(:balance, (nothing,), :eq,
            (emit, c, parameters, key) -> begin
                @coeff emit production[1] ← 1.0
                @coeff emit production[2] ← 1.0
            end,
            (c, parameters, key) -> c.demand[parameters.time_step]),),
        objective=Objective(:min, (emit, c, parameters) -> begin
            @coeff emit production[1] ← c.cost[parameters.time_step, 1]
            @coeff emit production[2] ← c.cost[parameters.time_step, 2]
        end),
    )
    closure_lp = LinearProblem(ConstData(closure_definition, const_data.data), p;
                               backend=HiGHSBackend(output_flag=false))
    @test closure_lp.const_data.data === lp.const_data.data
    @test isconcretetype(fieldtype(typeof(closure_lp), :numerical_updates))
end

@testset "Changing coefficient support" begin
    definition = ProblemDefinition(
        variables=(x=VariableBlock(ChangingSupportVariable),
                   y=VariableBlock(ChangingSupportVariable)),
        constraints=(cap=ConstraintBlock(:le, ChangingSupportCap),),
        objective=Objective(:max, ChangingSupportGain),
    )
    for backend in (HiGHSBackend(output_flag=false), ClpBackend())
        model = LinearProblem(ConstData(definition), ChangingSupportPoint(false); backend)
        @test objective_value(solve(model, ChangingSupportPoint(false))) ≈ 1.0
        @test objective_value(solve(model, ChangingSupportPoint(true))) ≈ 2.0
        @test objective_value(solve(model, ChangingSupportPoint(false))) ≈ 1.0
    end
end

@testset "Binary variables with a compatible backend" begin
    definition = ProblemDefinition(
        variables=(decision=VariableBlock(BinaryDecision; domain=:binary),),
        constraints=(floor=ConstraintBlock(:ge, BinaryFloor),),
        objective=Objective(:min, BinaryCost),
    )
    data = ConstData(definition)

    highs = LinearProblem(data, nothing; backend=HiGHSBackend(output_flag=false))
    @test solve(highs, nothing).variables.decision[nothing] == 1.0

    clp = LinearProblem(data, nothing; backend=ClpBackend())
    @test_throws ArgumentError solve(clp, nothing)
end

@testset "ConstData for transport and reservoir" begin
    arcs = [(1, 2), (2, 3), (1, 3)]
    tp = TransportParameters(1)
    transport_data = transport_const_data(
        arcs,
        1:3,
        [8.0 8.0 5.0],
        [1.0 1.0 3.0],
        [10.0 0.0 -10.0],
    )
    tlp = LinearProblem(transport_data, tp; backend=HiGHSBackend(output_flag=false))
    ts = solve(tlp, tp)
    @test fieldnames(TransportParameters) == (:pdt,)
    @test isbitstype(typeof(tp))
    @test transport_data.data.capacity isa Matrix{Float64}
    @test transport_data.data.injection isa Matrix{Float64}
    @test ts.variables.flow[(1, 2)] == 8.0
    @test ts.variables.flow[(2, 3)] == 8.0
    @test objective_value(ts) == 22.0

    rp = ReservoirParameters(1)
    reservoir_data = reservoir_const_data(
        [50.0],
        100.0,
        [20.0, 20.0, 20.0, 20.0],
        [100.0, 100.0, 100.0, 100.0],
        [1.0, 1.0, 2.0, 1.0],
        [5.0 10.0 4.0 2.0],
        [30.0 80.0 45.0 60.0],
        [0.9, 0.9, 0.9, 0.9],
        [50.0, 20.0],
        [0.0, 1500.0],
    )
    rlp = LinearProblem(reservoir_data, rp; backend=HiGHSBackend(output_flag=false))
    rs = solve(rlp, rp)
    @test fieldnames(ReservoirParameters) == (:scenario,)
    @test isbitstype(typeof(rp))
    @test reservoir_data.data.inflow isa Matrix{Float64}
    @test reservoir_data.data.price isa Matrix{Float64}
    @test termination_status(rs) == :optimal
    @test all(isfinite, rs.variables.stock.values)
    @test length(rs.variables.stock.keys) == 5
    @test rs.variables.stock[0] == 50.0
end

@testset "Clp warm coefficient update" begin
    data = dispatch_const_data(
        [100.0, 100.0],
        [80.0 60.0; 80.0 60.0],
        [40.0 70.0; 45.0 70.0],
    )
    p = DispatchParameters(1)
    changed_cost = DispatchParameters(2)
    lp = LinearProblem(data, Parameters(p); backend=ClpBackend())
    @test objective_value(solve(lp, p)) == 4600.0
    solution = solve(lp, changed_cost)
    @test solution.variables.production[1] == 80.0
    @test objective_value(solution) == 5000.0
end

@testset "Lazy solution materialization" begin
    lp = LinearProblem(dispatch_study(), DispatchParameters(1); backend=ClpBackend())

    first_solution = @inferred solve(lp, DispatchParameters(1))
    @test propertynames(first_solution) == (:variables, :constraints, :objective_value, :termination_status)
    @test getfield(first_solution, :variables_cache) === nothing
    @test getfield(first_solution, :constraints_cache) === nothing
    @test objective_value(first_solution) == 4600.0
    @test termination_status(first_solution) == :optimal
    @test getfield(first_solution, :variables_cache) === nothing
    @test getfield(first_solution, :constraints_cache) === nothing

    second_solution = solve(lp, DispatchParameters(2))
    @test second_solution === first_solution
    @test objective_value(first_solution) == 5300.0
    @test getfield(second_solution, :variables_cache) === nothing

    @test second_solution.variables.production[2] == 30.0
    saved_production = second_solution.variables.production
    @test getfield(second_solution, :variables_cache) !== nothing
    @test getfield(second_solution, :constraints_cache) === nothing

    third_solution = solve(lp, DispatchParameters(1))
    @test third_solution === second_solution
    @test saved_production[2] == 30.0
    @test getfield(third_solution, :variables_cache) === nothing
    @test third_solution.variables.production[2] == 20.0
    @test dual(third_solution, :balance) isa Float64
    @test getfield(third_solution, :constraints_cache) === nothing
    @test third_solution.constraints.balance[nothing] isa Float64
    @test getfield(third_solution, :variables_cache) !== nothing
    @test getfield(third_solution, :constraints_cache) !== nothing
    @test solve(lp, DispatchParameters(1)) === third_solution

    bellman_lp = LinearProblem(
        dispatch_study(), DispatchParameters(1); backend=ClpBackend(),
    )
    bellman_solution = solve(bellman_lp, DispatchParameters(1))
    @test @inferred(objective_value(bellman_lp, DispatchParameters(2))) == 5300.0
    @test bellman_lp.workspace.solution === bellman_solution
    @test objective_value(bellman_solution) == 5300.0
    @test getfield(bellman_solution, :variables_cache) === nothing
    same_point_solution = solve(bellman_lp, DispatchParameters(2))
    @test same_point_solution === bellman_solution
    @test objective_value(same_point_solution) == 5300.0
    @test solve(bellman_lp, DispatchParameters(2)) === same_point_solution
end

@testset "Selective dual access" begin
    for backend in (HiGHSBackend(output_flag=false), ClpBackend())
        lp = LinearProblem(dispatch_study(), DispatchParameters(1); backend)
        first_solution = solve(lp, DispatchParameters(1))
        first_dual = dual(first_solution, :balance)
        @test first_dual isa Float64
        @test getfield(first_solution, :constraints_cache) === nothing
        @test solve(lp, DispatchParameters(2)) === first_solution
        @test dual(first_solution, :balance) isa Float64
        @test getfield(first_solution, :constraints_cache) === nothing
    end
end

@testset "Finite and infinite bound transitions" begin
    struct BoundConstData
        lower::Vector{Float64}
        upper::Vector{Float64}
    end

    struct BoundParameters
        pdt::Int
    end

    definition = ProblemDefinition(
        variables=(x=VariableBlock(:x;
            lower=(c, p, _) -> c.lower[p.pdt],
            upper=(c, p, _) -> c.upper[p.pdt]),),
        constraints=(floor=ConstraintBlock(:floor, (nothing,), :ge,
            (emit, c, p, _) -> (@coeff emit x ← 1.0), (c, p, _) -> -10.0),),
        objective=Objective(:min, (emit, c, p) -> (@coeff emit x ← 1.0)),
    )
    data = ConstData(definition, BoundConstData([0.0, -Inf], [10.0, Inf]))
    lp = LinearProblem(data, BoundParameters(1);
                       backend=HiGHSBackend(output_flag=false))
    @test primal(solve(lp, BoundParameters(1)), :x) == 0.0
    @test solve(lp, BoundParameters(2)).variables.x[nothing] == -10.0
    x_index = lp.workspace.variable_indices[(:x, nothing)]
    @test lp.workspace.lower_values[x_index] == -Inf
    @test lp.workspace.upper_values[x_index] == Inf
    @test primal(solve(lp, BoundParameters(1)), :x) == 0.0
    @test lp.workspace.lower_values[x_index] == 0.0
    @test lp.workspace.upper_values[x_index] == 10.0
end

@testset "Finite and inactive constraint right-hand sides" begin
    struct RowParameters
        active::Bool
    end

    definition = ProblemDefinition(
        variables=(x=VariableBlock(:x; lower=0.0, upper=10.0),),
        constraints=(cap=ConstraintBlock(
            :cap,
            (nothing,),
            :le,
            (equation, data, parameters, _) -> (@coeff equation x ← 1.0),
            (data, parameters, _) -> parameters.active ? 3.0 : Inf,
        ),),
        objective=Objective(:max, (expression, data, parameters) ->
            (@coeff expression x ← 1.0)),
    )
    data = ConstData(definition)

    for backend in (HiGHSBackend(output_flag=false), ClpBackend())
        model = LinearProblem(data, RowParameters(true); backend)
        @test solve(model, RowParameters(true)).variables.x[nothing] == 3.0
        @test solve(model, RowParameters(false)).variables.x[nothing] == 10.0
        @test solve(model, RowParameters(true)).variables.x[nothing] == 3.0
        # Both direct backends use the GLP/HiGHS row-dual convention.
        @test dual(solve(model, RowParameters(true)), :cap) ≈ 1.0
    end
end

@testset "Batched RHS updates skip unchanged LHS" begin
    lhs_calls = Ref(0)
    definition = ProblemDefinition(
        variables=(x=VariableBlock(:x, 1:4; lower=0.0, upper=10.0),),
        constraints=(cap=ConstraintBlock(1:4, :le, BatchRow(lhs_calls)),),
        objective=Objective(:max, BatchRowObjective),
    )
    data = ConstData(definition)

    for backend in (HiGHSBackend(output_flag=false), ClpBackend())
        lhs_calls[] = 0
        model = LinearProblem(data, BatchRowPoint(0); backend)
        first_solution = solve(model, BatchRowPoint(0))
        @test lhs_calls[] == 4
        @test objective_value(first_solution) == 10.0

        second_solution = solve(model, BatchRowPoint(1))
        @test second_solution === first_solution
        @test lhs_calls[] == 4
        @test objective_value(second_solution) == 14.0
        @test all(second_solution.variables.x[key] == key + 1 for key in 1:4)

        @test objective_value(solve(model, BatchRowPoint(0))) == 10.0
        @test lhs_calls[] == 4
    end
end

@testset "Failed solve keeps the model and invalidates the cached solution" begin
    good = DispatchParameters(1)
    infeasible = DispatchParameters(4)
    for backend in (() -> HiGHSBackend(output_flag=false), ClpBackend)
        lp = LinearProblem(dispatch_study(), good; backend=backend())
        good_solution = solve(lp, good)
        @test objective_value(good_solution) == 4600.0
        solve_error = try
            solve(lp, infeasible)
            nothing
        catch error
            error
        end
        @test solve_error isa SolveError
        @test solve_error.status == :infeasible
        # GLP policy: the numerical model is kept, updated to the failed point.
        @test lp.workspace.built
        @test lp.workspace.parameter_snapshot == Parameters(infeasible)
        @test lp.workspace.solution === good_solution
        @test_throws ErrorException objective_value(good_solution)
        @test_throws ErrorException good_solution.variables.production[1]
        # Asking again for the failed point raises again, without a cached answer.
        @test_throws SolveError solve(lp, infeasible)
        @test_throws SolveError objective_value(lp, infeasible)
        # The next point is updated from the failed one.
        @test objective_value(solve(lp, good)) == 4600.0
        @test solve(lp, good) === good_solution
        @test objective_value(solve(lp, DispatchParameters(2))) ==
              objective_value(solve(LinearProblem(dispatch_study(), good; backend=backend()),
                                    DispatchParameters(2)))

        bellman_lp = LinearProblem(dispatch_study(), good; backend=backend())
        @test objective_value(bellman_lp, good) == 4600.0
        @test_throws SolveError objective_value(bellman_lp, infeasible)
        @test bellman_lp.workspace.solution === nothing
        @test bellman_lp.workspace.last_status === :infeasible
        @test objective_value(bellman_lp, good) == 4600.0
        @test bellman_lp.workspace.last_status === :optimal
    end

    # Infeasible from the very first solve: the model is built and kept too.
    lp = LinearProblem(dispatch_study(), infeasible; backend=HiGHSBackend(output_flag=false))
    @test_throws SolveError objective_value(lp, infeasible)
    @test lp.workspace.built
    @test objective_value(lp, good) == 4600.0
end

@testset "iteration_space and constraint_sense auto-wiring" begin
    struct AutoVar <: AbstractVariable end
    ReactiveModel.iteration_space(::AutoVar, data) = 1:data.n
    ReactiveModel.lower_bound(::AutoVar, data, params, i) = 0.0
    ReactiveModel.upper_bound(::AutoVar, data, params, i) = 10.0

    struct AutoConstraint <: AbstractConstraint end
    ReactiveModel.iteration_space(::AutoConstraint, data) = 1:data.n
    ReactiveModel.constraint_sense(::AutoConstraint) = :eq
    ReactiveModel.lhs_terms!(::AutoConstraint, emit, data, params, i) = @coeff emit x[i] ← 1.0
    ReactiveModel.rhs_value(::AutoConstraint, data, params, i) = Float64(i)

    struct AutoObj <: AbstractObjective end
    ReactiveModel.objective_terms!(::AutoObj, emit, data, params) = (@coeff emit x[1] ← 1.0)

    data = (n = 3,)
    vb = VariableBlock(:x, AutoVar(), data)
    cb = ConstraintBlock(:c, AutoConstraint(), data)
    @test vb.keys == [1, 2, 3]
    @test cb.keys == [1, 2, 3]
    @test cb.sense === :eq

    prob_def = ProblemDefinition(variables=(x=vb,), constraints=(c=cb,), objective=Objective(:min, AutoObj()))
    lp = LinearProblem(ConstData(prob_def, data), Parameters((s=1,)))
    sol = solve(lp, (s=1,))
    @test is_optimal(sol)
    @test primal(sol, :x, 1) == 1.0
    @test primal(sol, :x, 2) == 2.0
    @test primal(sol, :x, 3) == 3.0
end


@testset "@coeff macro and set_coeff! API" begin
    struct CoeffVar <: AbstractVariable end
    ReactiveModel.iteration_space(::CoeffVar, data) = 1:data.n
    ReactiveModel.variable_bounds(::CoeffVar) = (0.0, 10.0)

    struct CoeffConstraint <: AbstractConstraint end
    ReactiveModel.iteration_space(::CoeffConstraint, data) = 1:data.n
    ReactiveModel.constraint_sense(::CoeffConstraint) = :eq
    function ReactiveModel.lhs_terms!(::CoeffConstraint, eq, data, params, i)
        @coeff eq x[i] ← 1.0
    end
    ReactiveModel.rhs_value(::CoeffConstraint, data, params, i) = Float64(i * 2)

    struct CoeffBlockConstraint <: AbstractConstraint end
    ReactiveModel.iteration_space(::CoeffBlockConstraint, data) = [1]
    ReactiveModel.constraint_sense(::CoeffBlockConstraint) = :le
    function ReactiveModel.lhs_terms!(::CoeffBlockConstraint, eq, data, params, k)
        @coeff begin
            x[1] ← 1.0
            x[2] ← 1.0
        end
    end
    ReactiveModel.rhs_value(::CoeffBlockConstraint) = 10.0

    struct CoeffObj <: AbstractObjective end
    function ReactiveModel.objective_terms!(::CoeffObj, eq, data, params)
        @coeff x[1] ← 1.0
    end

    data = (n = 2,)
    vb = VariableBlock(:x, CoeffVar(), data)
    cb1 = ConstraintBlock(:c1, CoeffConstraint(), data)
    cb2 = ConstraintBlock(:c2, CoeffBlockConstraint(), data)

    prob_def = ProblemDefinition(
        variables=(x=vb,),
        constraints=(c1=cb1, c2=cb2),
        objective=Objective(:min, CoeffObj())
    )
    lp = LinearProblem(ConstData(prob_def, data), Parameters((s=1,)))
    sol = solve(lp, (s=1,))
    @test is_optimal(sol)
    @test primal(sol, :x, 1) == 2.0
    @test primal(sol, :x, 2) == 4.0
end


@testset "Declarative @model DSL" begin
    @model TestDSL(data, sel) begin
        R = 1:data.nR
        T = 1:data.nT

        @variable x[r in R, t in T] in [0.0, 10.0]
        @variable s                 in [0.0, Inf]

        @constraint Flow[r in R, t in T] begin
            x[r, t] <= data.cap[r, t]
        end

        @constraint Total begin
            sum(x[r, t] for r in R, t in T) + s == data.target[sel.step]
        end

        @objective Min begin
            100.0 * s + sum(data.cost[r] * x[r, t] for r in R, t in T)
        end
    end

    data = (
        nR = 2,
        nT = 2,
        cap = [5.0 5.0; 8.0 8.0],
        cost = [2.0, 5.0],
        target = [10.0, 15.0],
    )

    prob_def = TestDSL(data)
    @test haskey(prob_def.variables.by_name, :x)
    @test haskey(prob_def.variables.by_name, :s)
    @test haskey(prob_def.constraints.by_name, :Flow)
    @test haskey(prob_def.constraints.by_name, :Total)

    sel1 = (step = 1,)
    lp = LinearProblem(ConstData(prob_def, data), Parameters(sel1))
    sol1 = solve(lp, sel1)
    @test is_optimal(sol1)
    @test primal(sol1, :x, (1, 1)) == 5.0
    @test primal(sol1, :x, (1, 2)) == 5.0

    # Test automatically inferred change predicates
    total_ctr = prob_def.constraints.Total.constraint
    @test ReactiveModel.lhs_has_changed(total_ctr, (step=1,), (step=2,)) == false
    @test ReactiveModel.rhs_has_changed(total_ctr, (step=1,), (step=2,)) == true
    @test ReactiveModel.rhs_has_changed(total_ctr, (step=1,), (step=1,)) == false

    flow_ctr = prob_def.constraints.Flow.constraint
    @test ReactiveModel.lhs_has_changed(flow_ctr, (step=1,), (step=2,)) == false
    @test ReactiveModel.rhs_has_changed(flow_ctr, (step=1,), (step=2,)) == false

    # Test warm update and zero allocation
    sel2 = (step = 2,)
    sol2 = solve(lp, sel2)
    @test is_optimal(sol2)
    allocs = @allocated solve(lp, sel2)
    @test allocs == 0
end

@testset "Automatic AST dependency inference and param terminology" begin
    # 1. Test model defined with (data, param) without any manual update clauses
    @model AutoInferenceModel(data, param) begin
        @variable x in (0.0, data.cap[param.t])
        @variable y in (0.0, 10.0)

        @constraint balance x + 2.0 * y == data.demand[param.t, param.scenario]
        @constraint bound_y y <= 0.0

        @objective Min begin
            data.cost[param.t] * x + 3.0 * y
        end
    end

    study_data = (
        cap = [5.0, 8.0],
        demand = [4.0 6.0; 7.0 9.0],
        cost = [10.0, 20.0],
    )
    def = AutoInferenceModel(study_data)
    
    x_var = def.variables.x.variable
    y_var = def.variables.y.variable
    bal_ctr = def.constraints.balance.constraint
    by_ctr = def.constraints.bound_y.constraint
    obj = def.objective.objective

    # x depends on param.t
    @test ReactiveModel.has_changed(x_var, (t=1, scenario=1), (t=2, scenario=1)) == true
    @test ReactiveModel.has_changed(x_var, (t=1, scenario=1), (t=1, scenario=2)) == false
    # y is completely constant
    @test ReactiveModel.has_changed(y_var, (t=1, scenario=1), (t=2, scenario=2)) == false

    # balance LHS has no param -> false
    @test ReactiveModel.lhs_has_changed(bal_ctr, (t=1, scenario=1), (t=2, scenario=2)) == false
    # balance RHS depends on (param.t, param.scenario)
    @test ReactiveModel.rhs_has_changed(bal_ctr, (t=1, scenario=1), (t=2, scenario=1)) == true
    @test ReactiveModel.rhs_has_changed(bal_ctr, (t=1, scenario=1), (t=1, scenario=2)) == true
    @test ReactiveModel.rhs_has_changed(bal_ctr, (t=1, scenario=1), (t=1, scenario=1)) == false

    # bound_y is completely constant
    @test ReactiveModel.lhs_has_changed(by_ctr, (t=1, scenario=1), (t=2, scenario=2)) == false
    @test ReactiveModel.rhs_has_changed(by_ctr, (t=1, scenario=1), (t=2, scenario=2)) == false

    # objective depends on param.t
    @test ReactiveModel.has_changed(obj, (t=1, scenario=1), (t=2, scenario=1)) == true
    @test ReactiveModel.has_changed(obj, (t=1, scenario=1), (t=1, scenario=2)) == false

    # 2. Test symbolic clauses (:t, (:t, :scenario)) and lambda clause ((old, new) -> ...)
    @model SymbolicClauseModel(data, param) begin
        @variable x in (0.0, 10.0) update_when :t
        @constraint c1 begin
            x <= 5.0
        end update_rhs_when (:t, :scenario)
        @constraint c2 begin
            x >= 1.0
        end update_rhs_when ((old, new) -> old.t != new.t)
        @objective Min begin
            x
        end update_when false
    end

    def_sym = SymbolicClauseModel(study_data)
    @test ReactiveModel.has_changed(def_sym.variables.x.variable, (t=1, scenario=1), (t=2, scenario=1)) == true
    @test ReactiveModel.has_changed(def_sym.variables.x.variable, (t=1, scenario=1), (t=1, scenario=2)) == false

    @test ReactiveModel.rhs_has_changed(def_sym.constraints.c1.constraint, (t=1, scenario=1), (t=1, scenario=2)) == true
    @test ReactiveModel.rhs_has_changed(def_sym.constraints.c2.constraint, (t=1, scenario=1), (t=1, scenario=2)) == false
    @test ReactiveModel.rhs_has_changed(def_sym.constraints.c2.constraint, (t=1, scenario=1), (t=2, scenario=2)) == true

    @test ReactiveModel.has_changed(def_sym.objective.objective, (t=1, scenario=1), (t=2, scenario=2)) == false

    # Solve and verify zero allocation with auto-inferred model
    lp = LinearProblem(ConstData(def, study_data), (t=1, scenario=1))
    sol1 = solve(lp, (t=1, scenario=1))
    @test is_optimal(sol1)
    @test primal(sol1, :x) == 4.0

    sol2 = solve(lp, (t=2, scenario=1))
    @test is_optimal(sol2)
    @test primal(sol2, :x) == 7.0

    # A literal named tuple makes Julia 1.13 wrap this measurement in a fresh
    # closure, including its dispatch-cache setup. Measure the warmed solve
    # directly, with the same numerical transition and zero-byte assertion.
    allocation_point = (t=1, scenario=1)
    allocs = @allocated solve(lp, allocation_point)
    @test allocs == 0
end

@testset "rebuild! for mutable ConstData contract" begin
    mutable struct MutableStudy
        capacity::Vector{Float64}
        demand::Vector{Float64}
    end

    @model RebuildModel(data, param) begin
        G = 1:2
        @variable prod[g in G] in [0.0, data.capacity[g]]
        @constraint balance begin
            prod[1] + prod[2] == data.demand[param.step]
        end
        @objective Min begin
            10.0 * prod[1] + 20.0 * prod[2]
        end
    end

    study = MutableStudy([10.0, 10.0], [5.0, 15.0])
    def = RebuildModel(study)
    lp = LinearProblem(ConstData(def, study), Parameters((step=1,)))

    # 1. Initial solve at step 1 (demand = 5.0 -> prod[1] = 5.0, obj = 50.0)
    sol1 = solve(lp, (step=1,))
    @test is_optimal(sol1)
    @test objective_value(sol1) == 50.0
    @test primal(sol1, :prod, 1) == 5.0

    # Shared study arrays require a rebuild of each independent numerical model.
    other = similar(lp)
    @test objective_value(other, (step=1,)) == 50.0

    # 2. Mutate study data in-place: capacity of cheap generator drops to 2.0
    study.capacity[1] = 2.0
    study.demand[2] = 9.0

    # Calling rebuild!(lp) triggers full reconstruction honoring the contract
    rebuild!(lp)
    rebuild!(other)
    sol1_rebuilt = solve(lp, (step=1,))
    @test is_optimal(sol1_rebuilt)
    @test objective_value(sol1_rebuilt) == 80.0
    @test primal(sol1_rebuilt, :prod, 1) == 2.0
    @test primal(sol1_rebuilt, :prod, 2) == 3.0
    @test objective_value(other, (step=1,)) == 80.0
    # A different selector does not replace the rebuild required for changed data.
    @test objective_value(lp, (step=2,)) == 160.0
    @test objective_value(other, (step=2,)) == 160.0

    # 3. Test rebuild! with new raw data payload
    study2 = MutableStudy([8.0, 8.0], [4.0, 12.0])
    rebuild!(lp, study2)
    sol1_new = solve(lp, (step=1,))
    @test is_optimal(sol1_new)
    @test objective_value(sol1_new) == 40.0
    @test primal(sol1_new, :prod, 1) == 4.0

    # 4. Test rebuild! before any solve has occurred
    lp_fresh = LinearProblem(ConstData(def, study), Parameters((step=1,)))
    @test lp_fresh.built == false
    rebuild!(lp_fresh)
    @test lp_fresh.built == false
    sol_fresh = solve(lp_fresh, (step=1,))
    @test is_optimal(sol_fresh)
end

@testset "deepcopy prohibition and safe similar(lp) cloning" begin
    def = ProblemDefinition(
        variables=(x=VariableBlock(:x; lower=0.0, upper=10.0),),
        constraints=(cap=ConstraintBlock(
            :cap, (nothing,), :le,
            (eq, d, p, k) -> (@coeff eq x ← 1.0),
            (d, p, k) -> 7.0,
        ),),
        objective=Objective(:max, (eq, d, p) -> (@coeff eq x ← 1.0)),
    )
    lp = LinearProblem(ConstData(def), (step=1,); backend=HiGHSBackend(output_flag=false))
    solve(lp, (step=1,))

    # 1. Verify deepcopy is strictly prohibited and throws ArgumentError
    @test_throws ArgumentError deepcopy(lp)
    @test_throws ArgumentError deepcopy(lp.workspace)
    @test_throws ArgumentError deepcopy(lp.backend)

    # 2. Verify deepcopy of a user container wrapping a LinearProblem is also blocked
    struct UserProblemHolder
        id::Int
        problem::LinearProblem
    end
    holder = UserProblemHolder(42, lp)
    @test_throws ArgumentError deepcopy(holder)

    # 3. Verify similar(lp) creates an independent, isolated problem
    lp2 = similar(lp)
    @test lp2 isa LinearProblem
    @test lp2 !== lp
    @test lp2.const_data === lp.const_data
    @test lp2.backend !== lp.backend
    @test lp2.backend.inner !== lp.backend.inner # Distinct native C pointers!

    # 4. Both models can be solved independently at different parameters
    sol1 = solve(lp, (step=1,))
    sol2 = solve(lp2, (step=1,))
    @test is_optimal(sol1)
    @test is_optimal(sol2)
    @test primal(sol1, :x) == 7.0
    @test primal(sol2, :x) == 7.0

    # 5. copy(lp) forwards to similar(lp) safely
    lp_copy = copy(lp)
    @test lp_copy !== lp
    @test lp_copy.backend !== lp.backend
    @test lp_copy.backend.inner !== lp.backend.inner
end

@testset "Targeted setup_stmts filtering in @model DSL" begin
    @model FilteredTest(data, param) begin
        a = data.a
        b = data.b
        c = data.c
        @variable x in [0.0, a]
        @constraint c1 begin
            x <= b
        end
        @objective Min begin
            c * x
        end
    end

    data = (a = 10.0, b = 5.0, c = 2.0)
    def = FilteredTest(data)
    lp = LinearProblem(ConstData(def, data), (t=1,); backend=ClpBackend())
    sol = solve(lp, (t=1,))
    @test is_optimal(sol)
    @test objective_value(sol) == 0.0
    @test primal(sol, :x) == 0.0

    # Changing parameter
    sol2 = solve(lp, (t=2,))
    @test is_optimal(sol2)
end

@testset "ClpJLBackend: Extension Trigger" begin
    @test_throws ErrorException ClpJLBackend()
end

include("test_clp_backend.jl")
if Base.find_package("ClpJL") !== nothing
    include("test_clpjl.jl")
else
    @testset "Optional ClpJL integration (not installed)" begin
        @test_skip false
    end
end
include("test_highs_jll.jl")
include("test_dsl_domains.jl")
include("test_typed_plan.jl")
include("test_no_moi.jl")


include("test_publication_regressions.jl")
include("test_affine_dsl.jl")

include("test_optional_backends.jl")
include("test_backend_assembly.jl")
include("test_numeric_bounds.jl")
include("test_commercial_backend_contracts.jl")

include("test_reactivity_conditions.jl")

include("test_dsl_grammar.jl")
include("test_safe_updates.jl")
include("test_objective_composition.jl")
include("test_architecture.jl")

include("test_lifecycle.jl")
