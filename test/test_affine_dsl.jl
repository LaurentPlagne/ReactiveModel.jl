module AffineDSLTests
using Test
using ReactiveModel
import ReactiveModel as RM

@model Aliases(data, point) begin
    t = point.t
    selected = data.values[t]
    limit, demand = (selected + 1.0, selected)
    coefficient = selected / 2.0
    price = coefficient + 1.0
    @variable x in [0.0, limit]
    @constraint balance coefficient * x == demand
    @objective Min price * x
end

@model TwoSided(data, point) begin
    offset = data.offset[point.t]
    @variable x in [0.0, 20.0]
    @variable y in [1.0, 1.0]
    @constraint balance 2.0 * (x + offset) == x / 2.0 + y + 5.0
    @objective Max -(x + y) / 2.0 + 2.0 * (offset + 1.0)
end

@model ConstantOnlyUpdate(data, point) begin
    bonus = data.bonus[point.t]
    @variable x in [3.0, 3.0]
    @objective Min x + bonus
end

@model SummedConstants(data, point) begin
    @variable x in [0.0, 20.0]
    @constraint balance sum(x / 2.0 + point.shift for i in 1:data.n if i > data.skip) == point.rhs
    @objective Min sum(point.choose ? (x + point.shift) : (2.0 * x - point.shift) for i in 1:data.n if i > data.skip) + 7.0
end

@model RightHandVariables(data, point) begin
    @variable x in [0.0, 20.0]
    @constraint cap point.offset + 8.0 >= point.a * x + 2.0
    @constraint floor point.a * x - point.offset >= 0.0
    @objective Max x + point.offset
end

struct Collision
    id::Int
end
Base.hash(::Collision, seed::UInt) = hash(0, seed)

module OwnerA
using ReactiveModel
@model Shared(data, point) begin
    @variable x in [0.0, 10.0]
    @constraint cap x <= 1.0
    @objective Max x + 10.0
end
end
module OwnerB
using ReactiveModel
@model Shared(data, point) begin
    @variable x in [0.0, 20.0]
    @constraint cap x <= 9.0
    @objective Min -x + 20.0
end
end

problem(model, data, point; backend=HiGHSBackend()) =
    LinearProblem(ConstData(model(data), data), point; backend)

@testset "Transitive dependencies across every numeric component" begin
    data = (values=[2.0, 6.0],)
    lp = problem(Aliases, data, (t=1, irrelevant=0))
    def = lp.definition
    p, q = (t=1, irrelevant=0), (t=2, irrelevant=0)
    for object in (def.variables.x.variable, def.objective.objective)
        @test has_changed(object, p, q)
        @test !has_changed(object, p, (t=1, irrelevant=1))
    end
    @test lhs_has_changed(def.constraints.balance.constraint, p, q)
    @test rhs_has_changed(def.constraints.balance.constraint, p, q)
    for t in (1, 2, 1)
        point = (t=t, irrelevant=0)
        sol = solve(lp, point)
        @test objective_value(sol) ≈ data.values[t] + 2.0
        @test primal(sol, :x) ≈ 2.0
        @test objective_value(sol) ≈ objective_value(problem(Aliases, data, point), point)
        @test lp.upper_values[1] == data.values[t] + 1.0
    end
    @test RM._formula_param_deps(:alias, [:(alias = identity(point))], :point) === :opaque
end

@testset "Affine normalization, updates and result paths" begin
    for backend in (HiGHSBackend, ClpBackend)
        data = (offset=[0.0, 1.5],)
        lp = problem(TwoSided, data, (t=1,); backend=backend())
        for t in (1, 2, 1)
            sol = solve(lp, (t=t,))
            x = (6.0 - 2.0 * data.offset[t]) / 1.5
            expected = -(x + 1.0) / 2.0 + 2.0 * (data.offset[t] + 1.0)
            @test primal(sol, :x) ≈ x
            @test objective_value(sol) ≈ expected
            @test objective_value(lp, (t=t,)) ≈ expected
            @test objective_value(sol) ≈ objective_value(problem(TwoSided, data, (t=t,); backend=backend()), (t=t,))
        end

        data = (bonus=[5.0, -2.0],)
        lp = problem(ConstantOnlyUpdate, data, (t=1,); backend=backend())
        sol = solve(lp, (t=1,))
        @test objective_value(sol) == 8.0
        @test solve(lp, (t=2,)) === sol
        @test objective_value(sol) == 1.0
        @test objective_value(lp, (t=1,)) == 8.0
        @test objective_value(sol) == 8.0
        @test objective_value(similar(lp), (t=2,)) == 1.0
        data.bonus[1] = 11.0
        rebuild!(lp)
        @test_throws ErrorException objective_value(sol)
        @test objective_value(solve(lp, (t=1,))) == 14.0
    end

    data = (n=3, skip=1)
    p = (shift=1.0, rhs=6.0, choose=true)
    lp = problem(SummedConstants, data, p)
    for point in (p, (shift=2.0, rhs=10.0, choose=false), p)
        sol = solve(lp, point)
        x = point.rhs - 2.0 * point.shift
        expected = 2.0 * (point.choose ? x + point.shift : 2.0 * x - point.shift) + 7.0
        @test primal(sol, :x) ≈ x
        @test objective_value(sol) ≈ expected
        @test objective_value(sol) ≈ objective_value(problem(SummedConstants, data, point), point)
    end
    empty_data = (n=3, skip=3)
    empty_point = (shift=1.0, rhs=0.0, choose=true)
    @test objective_value(problem(SummedConstants, empty_data, empty_point), empty_point) == 7.0

    lp = problem(RightHandVariables, nothing, (a=1.0, offset=0.0))
    for p in ((a=1.0, offset=0.0), (a=2.0, offset=2.0), (a=1.0, offset=0.0))
        sol = solve(lp, p)
        @test primal(sol, :x) ≈ (6.0 + p.offset) / p.a
        @test objective_value(sol) ≈ (6.0 + p.offset) / p.a + p.offset
    end
end

@testset "Colliding keys retain identity through updates and collections" begin
    keys = [Collision(1), Collision(2)]
    for keyset in (keys, UInt[1, 2], [(:a, 1), (:a, 2)])
        def = ProblemDefinition(
            variables=(x=VariableBlock(:x, keyset; lower=0.0, upper=10.0),),
            constraints=(fix=ConstraintBlock(:fix, keyset, :eq,
                (eq, d, p, k) -> set_coeff!(eq, :x, k, 1.0),
                (d, p, k) -> findfirst(isequal(k), keyset) * p.scale),),
            objective=Objective(:min, (eq, d, p) -> begin
                for key in keyset
                    set_coeff!(eq, :x, key, 1.0)
                end
            end))
        lp = LinearProblem(def, (scale=1.0,))
        for scale in (1.0, 2.0, 1.0)
            sol = solve(lp, (scale=scale,))
            @test objective_value(sol) == 3.0 * scale
            @test primal(sol, :x, keyset[1]) == scale
            @test sol.variables.x[keyset[2]] == 2.0 * scale
            @test dual(sol, :fix, keyset[1]) == 1.0
            @test activity(sol, :fix, keyset[2]) == 2.0 * scale
        end
    end
    a, b = keys
    def = ProblemDefinition(variables=(x=VariableBlock(:x, keys; lower=0.0, upper=1.0),),
                            constraints=NamedTuple(),
                            objective=Objective(:min, (eq, d, p) -> set_coeff!(eq, :x, a, Inf)))
    @test_throws ArgumentError objective_value(LinearProblem(def, nothing), nothing)
    @test hash(RM.VarKey(:x, a)) == hash(RM.VarKey(:x, b))
    @test length(Dict(RM.VarKey(:x, a)=>1, RM.VarKey(:x, b)=>2)) == 2
    @test objective_value(problem(OwnerA.Shared, nothing, nothing), nothing) == 11.0
    @test objective_value(problem(OwnerB.Shared, nothing, nothing), nothing) == 11.0
    da, db = OwnerA.Shared(nothing), OwnerB.Shared(nothing)
    @test typeof(da.constraints.cap.constraint) !== typeof(db.constraints.cap.constraint)
    @test typeof(da.objective.objective) !== typeof(db.objective.objective)
end

@testset "Unsupported mathematics fails explicitly" begin
    vars = Set([:x, :y])
    for expr in (:(x * y), :(x / y), :(sin(x)), :(x^2), :(x[y]), :(x > 0 ? x : y))
        @test_throws ArgumentError RM.lower_affine_expr(expr, vars)
    end
    for expr in (:(x < 2), :(x > 2))
        @test_throws ArgumentError RM.split_relation(expr)
        @test_throws ArgumentError RM.parse_variable_bounds(expr)
    end
    @test_throws ArgumentError RM._parse_when_args([:update_when])
    @test_throws ArgumentError RM._parse_when_args([:unknown, true])
    for body in (quote
            @variable x
            @constraint c begin x == 1; x == 2 end
        end, quote
            @variable x
            @objective Unknown x
        end, quote
            @variable x
            @objective Min begin x; x + 1 end
        end, quote
            @variable x
            @objective Min x
            @objective Max x
        end)
        @test_throws ArgumentError RM._build_model_expr(:Invalid, :data, :point, body)
    end
end
end
