# Regression checks included by test/runtests.jl.
module ReactiveModelPublicationAudit

using Test
using ReactiveModel
import ReactiveModel as RM

@model AliasedDemand(data, point) begin
    demand = data.demand[point.t]
    @variable x in [0.0, 10.0]
    @constraint balance x == demand
    @objective Min x
end

@model AffineConstant(data, point) begin
    @variable x in [0.0, 10.0]
    @constraint balance x + 1.0 == 3.0
    @objective Min x
end

@model ObjectiveConstant(data, point) begin
    @variable x in [0.0, 10.0]
    @constraint balance x == 3.0
    @objective Min x + 5.0
end

@model BinaryDomain(data, point) begin
    @variable x in [0.0, 1.0] domain = :binary
    @constraint floor x >= 0.5
    @objective Min x
end

@model CachedScalarDual(data, point) begin
    @variable x in [0.0, 10.0]
    @constraint floor x >= 2.0
    @objective Min x
end

@model RuntimeShape(data, point) begin
    @variable x[i in 1:data.count] in [0.0, 10.0]
    @constraint floor x[1] >= 2.0
    @objective Min x[1]
end

struct CollidingKey
    id::Int
end
Base.hash(::CollidingKey, seed::UInt) = hash(0, seed)

module OwnerA
using ReactiveModel
@model SameModel(data, point) begin
    @variable x in [0.0, 1.0]
    @objective Min x
end
end

module OwnerB
using ReactiveModel
@model SameModel(data, point) begin
    @variable x in [0.0, 9.0]
    @objective Min x
end
end

make_problem(definition, data=nothing, point=nothing) =
    LinearProblem(ConstData(definition, data), point; backend=HiGHSBackend())

function run()
    @testset "Audit : régressions à corriger avant publication" begin
        @testset "Dépendance indirecte du RHS" begin
            data = (demand=[2.0, 5.0],)
            definition = AliasedDemand(data)
            lp = make_problem(definition, data, (t=1,))
            @test objective_value(lp, (t=1,)) ≈ 2.0
            warm = objective_value(lp, (t=2,))
            cold = objective_value(make_problem(definition, data, (t=2,)), (t=2,))
            println("Alias RHS : warm=", warm, ", cold=", cold)
            @test cold ≈ 5.0
            @test warm ≈ cold
        end

        @testset "Constantes affines" begin
            lhs = objective_value(make_problem(AffineConstant(nothing)), nothing)
            obj = objective_value(make_problem(ObjectiveConstant(nothing)), nothing)
            println("Constante LHS : obtenu=", lhs, ", attendu=2.0")
            println("Constante objectif : obtenu=", obj, ", attendu=8.0")
            @test lhs ≈ 2.0
            @test obj ≈ 8.0
        end

        @testset "Domaine binaire du DSL" begin
            definition = BinaryDomain(nothing)
            value = objective_value(make_problem(definition), nothing)
            println("Domaine DSL : ", definition.variables.x.domain, ", optimum=", value)
            @test definition.variables.x.domain === :binary
            @test value ≈ 1.0
        end

        @testset "Dual scalaire indépendant du cache" begin
            sol = solve(make_problem(CachedScalarDual(nothing)), nothing)
            before = dual(sol, :floor)
            sol.constraints
            after = try
                dual(sol, :floor)
            catch err
                err
            end
            println("Dual scalaire : avant=", before, ", après=", after)
            @test before ≈ 1.0
            @test isequal(before, after)
        end

        @testset "Identité des clés et des modèles" begin
            a, b = CollidingKey(1), CollidingKey(2)
            @test a != b
            @test hash(a) == hash(b)
            @test RM.VarKey(:x, a) != RM.VarKey(:x, b)
            keys = [a, b]
            definition = ProblemDefinition(
                variables=(x=VariableBlock(:x, keys;
                    lower=(d, p, k) -> Float64(k.id),
                    upper=(d, p, k) -> Float64(k.id)),),
                constraints=NamedTuple(),
                objective=Objective(:min, (emit, d, p) -> begin
                    for k in keys
                        set_coeff!(emit, :x, k, 1.0)
                    end
                end),
            )
            collision_value = objective_value(make_problem(definition), nothing)
            println("Clés distinctes de même hash : optimum=", collision_value, ", attendu=3.0")
            @test collision_value ≈ 3.0
            va = OwnerA.SameModel(nothing).variables.x.variable
            vb = OwnerB.SameModel(nothing).variables.x.variable
            println("Modèles de modules distincts : même type=", typeof(va) === typeof(vb),
                    ", bornes du modèle A=", variable_bounds(va, nothing, nothing, 0))
            @test typeof(va) !== typeof(vb)
            @test variable_bounds(va, nothing, nothing, 0) == (0.0, 1.0)
        end

        @testset "Clonage avant le premier solve" begin
            data = (demand=[2.0, 5.0],)
            lp = make_problem(AliasedDemand(data), data, (t=1,))
            clone = similar(lp)
            println("Type paramètres : original=", typeof(lp).parameters[2],
                    ", clone=", typeof(clone).parameters[2])
            @test typeof(clone) === typeof(lp)
        end

        @testset "rebuild! conserve les espaces d'indices (contrat numérique)" begin
            data = (count=1,)
            lp = make_problem(RuntimeShape(data), data)
            solve(lp, nothing)
            rebuild!(lp, (count=2,))
            println("rebuild! : données count=", lp.const_data.data.count,
                    ", nombre de clés=", length(variable_block_keys(lp, :x)))
            @test length(variable_block_keys(lp, :x)) == 1
        end
    end
end

end # module

ReactiveModelPublicationAudit.run()
