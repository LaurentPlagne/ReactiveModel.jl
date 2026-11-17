module DSLGrammarTests
using Test
using ReactiveModel
const RM = ReactiveModel

@model Grammar(data, point) begin
    indices = data.indices
    rate, shift = (point.rate, point.shift)
    @variable x[i ∈ indices; i > 0] in [0.0, 20.0]
    @variable y[i = indices, j in 1:i; i > 0 && j > 0] in [0.0, 1.0]
    @variable anchor >= 0.0
    @constraint fixed[i in indices; i > 0] rate * (x[i] + shift) == point.rhs + x[i] / 2.0
    @constraint zero anchor == 0.0
    @objective Min sum(point.choose ? x[i] + shift : 2.0 * x[i] - shift for i in indices if i > 0) +
        sum(sum(y[i,j] for j in 1:i) for i in indices if i > 0) + 3.0
end

# Argument and iterator names must not capture compiler temporaries.
@model Names(key, eq) begin
    @variable x in [0.0, eq.limit]
    @variable y[emit in 1:key.n] in [0.0, 1.0]
    @constraint floor x >= eq.rhs
    @objective Min x + sum(y[emit] for emit in 1:key.n)
end

bounds_pair(data, point) = (data.lower, point.cap)
@model BoundForms(data, point) begin
    @variable free
    @variable upper <= point.cap
    @variable lower >= point.floor
    @variable 0.0 <= ranged <= point.cap
    @variable pair in (data.lower, point.cap)
    @variable helper in bounds_pair(data, point)
    @constraint fix free == point.floor
    @objective Min free - upper + lower + ranged + pair + helper
end

@model ExternalRows(data, point) begin
    @external_variables x
    @constraint cap x <= point.cap
    @objective Max x
end

@testset "Documented DSL grammar: numeric and structural behavior" begin
    for backend in (HiGHSBackend, ClpBackend), indices in ([0, 1, 2], Int[])
        data = (indices=indices,)
        p = (rate=1.0, shift=1.0, rhs=3.0, choose=true)
        lp = LinearProblem(ConstData(Grammar(data), data), p; backend=backend())
        for q in (p, (rate=2.0, shift=0.5, rhs=4.0, choose=false), p)
            sol = solve(lp, q)
            x = (q.rhs - q.rate*q.shift) / (q.rate - 0.5)
            count = length(filter(>(0), indices))
            expected = count*(q.choose ? x + q.shift : 2x - q.shift) + 3.0
            @test objective_value(sol) ≈ expected
            @test length(lp.definition.variables.x.keys) == count
            @test length(lp.definition.variables.y.keys) == sum((i for i in indices if i > 0); init=0)
            fresh = LinearProblem(ConstData(Grammar(data), data), q; backend=backend())
            @test objective_value(sol) ≈ objective_value(fresh, q)
            for i in indices
                i > 0 && (@test primal(sol, :x, i) ≈ x)
            end
        end
    end
    data = (n=3,)
    p = (limit=10.0, rhs=2.0)
    lp = LinearProblem(ConstData(Names(data), data), p)
    @test objective_value(lp, p) == 2.0
    @test objective_value(lp, (limit=9.0, rhs=4.0)) == 4.0
end

@testset "Bounds and external-family composition" begin
    data = (lower=1.0,)
    for backend in (HiGHSBackend, ClpBackend)
        p = (cap=10.0, floor=2.0)
        lp = LinearProblem(ConstData(BoundForms(data), data), p; backend=backend())
        for q in (p, (cap=8.0, floor=3.0), p)
            @test objective_value(lp, q) ≈ 2*q.floor-q.cap+2*data.lower
        end
        rows = ExternalRows(nothing)
        definition = ProblemDefinition(variables=(x=VariableBlock(:x; lower=0.0,upper=10.0),),
            constraints=rows.constraints, objective=rows.objective)
        composed = LinearProblem(definition, (cap=4.0,); backend=backend())
        @test objective_value(composed, (cap=4.0,)) == 4.0
        @test objective_value(composed, (cap=7.0,)) == 7.0
    end
end

@testset "Unsupported DSL grammar fails during expansion" begin
    invalid = [
        "@variable x[I]", "@variable x[]", "@variable x[(i,j) in data.pairs]",
        "@variable x[i in 1:2, i in 1:2]", "@variable x in [0, 1, 2]",
        "@variable x in []", "@variable x[i in 1:2; i > 0; i < 2]",
        "@variable x[i in 1:point.n]", "@variable x[i in 1:2; i > point.skip]",
        "n = point.n; @variable x[i in 1:n]", "p = identity(point); @variable x[i in 1:p.n]",
        "@variable x domain=point.domain", "@variable x; @variable y in [0, x]",
        "@variable x; alias = x + 1; @objective Min alias",
        "@variable x; @variable x", "@constraint c 1 == 1; @constraint c 2 == 2",
        "@variable", "@constraint c", "@objective Min", "@external_variables",
        "@external_variables (x, 1)", "@external_variables ()", "a::Int = 1", "@unknown x", "@variable x; @constraint c (x = 1)",
        "@variable x; @constraint c 0 <= x <= 1", "@variable x; @objective Min x^2",
        "@variable x; @objective Min sum(x for i in 1:2 for j in 1:2)",
        "for i in 1:2; @variable x; end", "if true; @variable x; end",
        "data.n = 2", "data[1] = 2", "f(x) = x + 1", "a = 1; a = 2",
        "a = b; b = 1", "a = a + 1", "point = 1", "a = begin 1; 2 end",
        "@variable x; @constraint c (x < 1)",
    ]
    for body in invalid
        expr = Meta.parse("begin $body end")
        @test_throws ArgumentError RM._build_model_expr(gensym(:Invalid), :data, :point, expr)
    end
    for signature in (:(Bad(a,b,c)), :(Bad(a,a)), :(Bad(a::Int,b)))
        expression = Expr(:macrocall, Symbol("@model"), LineNumberNode(1), signature,
            quote @variable x end)
        try
            macroexpand(@__MODULE__, expression)
            @test false
        catch error
            @test error isa ArgumentError || (error isa LoadError && error.error isa ArgumentError)
        end
    end
end
end
