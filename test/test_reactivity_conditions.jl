using Test
using ReactiveModel

module ReactivityConditionsTests
using Test
using ReactiveModel
const RM = ReactiveModel

@model Explicit(data, point) begin
    @variable x in [0.0, point.limit] update_when :limit
    @constraint floor point.a * (x + point.shift) >= point.rhs update_lhs_when :a update_rhs_when (:a, :shift, :rhs)
    @objective Min point.cost * x + point.offset update_when ((previous, current) -> !isequal((previous.cost, previous.offset), (current.cost, current.offset)))
end

@model SignedZero(data, point) begin
    @variable x in [0.0, signbit(point.zero) ? 2.0 : 1.0]
    @objective Max x
end

@model Compatibility(data, point) begin
    alias = point.t
    @variable x in [0.0, 10.0] when :t
    @constraint lambda x <= point.t update_rhs_when ((n, o) -> !isequal(n.t, o.t))
    @constraint legacy x <= point.t rhs_when (o.t != n.t)
    @constraint shorthand x <= alias update_rhs_when changed(alias)
    @constraint general_first x <= point.t update_when false update_rhs_when :t
    @constraint specific_first x <= point.t update_rhs_when :t update_when false
    @objective Max x changed_when false
end

@model Opaque(data, point) begin
    @variable x in [0.0, identity(point).limit]
    @objective Max x
end

@testset "Explicit reactivity: incremental equals fresh" begin
    data = nothing
    points = (
        (limit=10.0, a=1.0, shift=1.0, rhs=3.0, cost=2.0, offset=7.0, irrelevant=0),
        (limit=9.0, a=2.0, shift=0.5, rhs=7.0, cost=4.0, offset=9.0, irrelevant=0),
        (limit=8.0, a=2.0, shift=1.0, rhs=6.0, cost=3.0, offset=2.0, irrelevant=1),
        (limit=8.0, a=2.0, shift=1.0, rhs=6.0, cost=3.0, offset=5.0, irrelevant=1),
        (limit=8.0, a=2.0, shift=1.0, rhs=6.0, cost=3.0, offset=5.0, irrelevant=2),
    )
    definition = Explicit(data)
    for backend in (HiGHSBackend, ClpBackend)
        lp = LinearProblem(ConstData(definition, data), first(points); backend=backend())
        for p in (points..., first(points))
            sol = solve(lp, p)
            fresh = LinearProblem(ConstData(definition, data), p; backend=backend())
            expected = p.rhs / p.a - p.shift
            @test isapprox(primal(sol, :x), expected)
            @test isapprox(objective_value(sol), p.cost * expected + p.offset)
            @test isapprox(objective_value(sol), objective_value(fresh, p))
        end
    end
    p, q = points[4:5]
    @test !has_changed(definition.variables.x.variable, p, q)
    @test !lhs_has_changed(definition.constraints.floor.constraint, p, q)
    @test !rhs_has_changed(definition.constraints.floor.constraint, p, q)
    @test !has_changed(definition.objective.objective, p, q)
end

@testset "Comparison contract and compatibility" begin
    def = Compatibility(nothing)
    p, q = (t=1, ignored=0), (t=2, ignored=0)
    @test has_changed(def.variables.x.variable, p, q)
    for name in (:lambda, :legacy, :shorthand, :general_first, :specific_first)
        tag = getproperty(def.constraints, name).constraint
        @test rhs_has_changed(tag, p, q)
        @test !rhs_has_changed(tag, p, (t=1, ignored=1))
        @test !lhs_has_changed(tag, p, q)
    end
    @test !has_changed(def.objective.objective, p, q)
    lp = LinearProblem(ConstData(def), p)
    @test objective_value(lp, p) == 1.0
    @test objective_value(lp, q) == 2.0
    @test objective_value(lp, p) == 1.0

    signed = SignedZero(nothing)
    tag = signed.variables.x.variable
    @test has_changed(tag, (zero=0.0,), (zero=-0.0,))
    @test !has_changed(tag, (zero=NaN,), (zero=NaN,))
    for backend in (HiGHSBackend, ClpBackend)
        lp = LinearProblem(ConstData(signed), (zero=0.0,); backend=backend())
        @test objective_value(lp, (zero=0.0,)) == 1.0
        @test objective_value(lp, (zero=-0.0,)) == 2.0
        @test objective_value(lp, (zero=0.0,)) == 1.0
    end
    opaque = Opaque(nothing).variables.x.variable
    @test has_changed(opaque, (limit=0.0,), (limit=-0.0,))
    @test !has_changed(opaque, (limit=NaN,), (limit=NaN,))
end

@testset "Malformed conditions are rejected at macro expansion" begin
    function model_with(declaration)
        RM._build_model_expr(gensym(:InvalidCondition), :data, :point,
            Expr(:block, declaration))
    end
    for condition in (:(()), :((:t, 2)), :((:t, point.rhs)),
                      :(previous -> true), :((previous, current, extra) -> true),
                      :((previous, previous) -> true), :(changed()), :(δ(point.t, point.rhs)),
                      1, nothing)
        declaration = Expr(:macrocall, Symbol("@variable"), LineNumberNode(1),
            :(x in [0.0, 10.0]), :update_when, condition)
        @test_throws ArgumentError model_with(declaration)
    end
    for clause in (:update_lhs_when, :update_rhs_when, :lhs_when, :rhs_when)
        @test_throws ArgumentError model_with(Expr(:macrocall, Symbol("@variable"),
            LineNumberNode(1), :(x in [0.0, 10.0]), clause, false))
        @test_throws ArgumentError model_with(Expr(:macrocall, Symbol("@objective"),
            LineNumberNode(1), :Min, 0.0, clause, false))
    end
    for clauses in ((:update_when, :t, :when, :rhs),
                    (:update_rhs_when, :t, :rhs_when, :rhs),
                    (:update_lhs_when, :t, :lhs_when, false),
                    (:update_when,))
        @test_throws ArgumentError model_with(Expr(:macrocall, Symbol("@constraint"),
            LineNumberNode(1), :cap, :(1.0 <= point.rhs), clauses...))
    end
end
end # module
