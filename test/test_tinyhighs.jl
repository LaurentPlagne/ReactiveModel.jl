module TinyHiGHSExtensionTests
using Test
using ReactiveModel
import TinyHiGHS

struct CollidingKey
    id::Int
end
Base.hash(::CollidingKey, seed::UInt) = hash(0, seed)

@model OptionalTinyModel(data, point) begin
    capacity = data.capacity[point.t]
    @variable x[k in data.keys] in [0.0, capacity]
    @constraint floor[k in data.keys] point.a * (x[k] + point.shift) >= point.rhs + k.id
    @objective Min sum(point.cost * x[k] + point.shift for k in data.keys)
end

@model OptionalTinySupport(data, point) begin
    @variable x in [0.0, 1.0]
    @variable y in [0.0, 1.0]
    @constraint cap x + (point.extra ? y : 0.0) <= 1.0
    @objective Max x + (point.extra ? 2.0 * y : 0.0) + point.bonus
end

include("trim/workload.jl")

@testset "TinyHiGHS optional extension" begin
    @test Base.get_extension(ReactiveModel, :ReactiveModelTinyHiGHSExt) !== nothing
    data = (capacity=[20.0, 30.0], keys=[CollidingKey(1), CollidingKey(2)])
    definition = ConstData(OptionalTinyModel(data), data)
    first_point = (t=1, a=1.0, shift=0.0, rhs=1.0, cost=2.0)
    next_point = (t=2, a=2.0, shift=1.0, rhs=4.0, cost=3.0)
    lp = LinearProblem(definition, first_point; backend=TinyHiGHSBackend())
    for point in (first_point, next_point, first_point)
        solution = @inferred solve(lp, point)
        reference = solve(LinearProblem(definition, point; backend=HiGHSBackend()), point)
        @test objective_value(solution) ≈ objective_value(reference)
        for key in data.keys
            @test primal(solution, :x, key) ≈ primal(reference, :x, key)
            @test dual(solution, :floor, key) ≈ dual(reference, :floor, key)
        end
        @test solution.variables.x[data.keys[1]] ≈ primal(reference, :x, data.keys[1])
    end
    clone = similar(lp)
    @test typeof(clone) === typeof(lp)
    @test clone.backend !== lp.backend
    @test objective_value(clone, next_point) ≈ objective_value(lp, next_point)
    rebuild!(lp)
    @test objective_value(lp, first_point) ≈ 10.0
    mktempdir() do dir
        file = joinpath(dir, "tiny.lp")
        backend_write_to_file(lp.backend, file)
        @test isfile(file)
        @test filesize(file) > 0
    end
    support = LinearProblem(ConstData(OptionalTinySupport(nothing)),
                            (extra=false, bonus=0.0); backend=TinyHiGHSBackend())
    for (point, expected) in (((extra=false, bonus=0.0), 1.0),
                              ((extra=true, bonus=0.0), 2.0),
                              ((extra=true, bonus=5.0), 7.0),
                              ((extra=false, bonus=0.0), 1.0))
        @test objective_value(support, point) == expected
    end
    @test ReactiveModelTrimWorkload.verify_lp(1; backend=TinyHiGHSBackend())
    @test ReactiveModelTrimWorkload.verify_lp(5; backend=TinyHiGHSBackend())
    binary = ProblemDefinition(
        variables=(x=VariableBlock(:x; lower=0.0, upper=1.0, domain=:binary),),
        constraints=NamedTuple(),
        objective=Objective(:min, (eq, d, p) -> set_coeff!(eq, :x, 1.0)),
    )
    @test_throws ArgumentError objective_value(LinearProblem(binary, nothing; backend=TinyHiGHSBackend()), nothing)
end
end
