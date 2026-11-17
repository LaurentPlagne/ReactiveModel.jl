module NumericBoundsTests

using Test, ReactiveModel
const RM = ReactiveModel

pair_bounds(data, point) = (point.lo, point.hi)

@model HelperBounds(data, point) begin
    @variable x[i in 1:data.n] in pair_bounds(data, point)
    @constraint floor[i in 1:data.n] x[i] >= -100.0
    @objective Max sum(point.cost * x[i] for i in 1:data.n)
end

@model LiteralBounds(data, point) begin
    @variable x[i in 1:data.n] in [point.lo, point.hi]
    @constraint floor[i in 1:data.n] x[i] >= -100.0
    @objective Max sum(point.cost * x[i] for i in 1:data.n)
end

struct FullBounds <: AbstractVariable end
RM.variable_bounds(::FullBounds, data, point, key) = (point.lo, point.hi)
struct CompactBounds{T} <: AbstractVariable
    lo::T
    hi::T
end
RM.variable_bounds(v::CompactBounds) = (v.lo, v.hi)
struct KeyBounds{T} <: AbstractVariable end
RM.variable_bounds(::KeyBounds{T}, key) where {T} = (T(key - 1), T(key + 4))
struct SeparateBounds <: AbstractVariable end
RM.lower_bound(::SeparateBounds, data, point, key) = point.lo
RM.upper_bound(::SeparateBounds, data, point, key) = point.hi
struct FaultyBounds <: AbstractVariable end
RM.variable_bounds(::FaultyBounds, data, point, key) =
    point.bad ? (nothing,point.hi) : (point.lo,point.hi)

function functional(variable, n; domain=:continuous)
    definition = ProblemDefinition(
        variables=(x=VariableBlock(:x, 1:n, variable; domain),),
        constraints=(floor=ConstraintBlock(:floor, 1:n, :ge,
            (e,d,p,k)->set_coeff!(e,:x,k,1.0), (d,p,k)->-100.0),),
        objective=Objective(:max, (e,d,p)->foreach(k->set_coeff!(e,:x,k,p.cost), 1:d.n)),
    )
    return ConstData(definition, (n=n,))
end

function sweep!(lp, points)
    for point in points
        objective_value(lp, point)
    end
    return nothing
end
function allocated_sweep(lp, points)
    sweep!(lp, points)
    return @allocated sweep!(lp, points)
end

@testset "Numeric bound conversion is inferred and precedes domain policy" begin
    for (lo, hi) in ((0,10), (0.0f0,10.0f0), (0,10.0f0), (1//3,7//2), (-Inf,Inf))
        @test @inferred(RM._domain_bounds(:continuous,lo,hi)) === (Float64(lo),Float64(hi))
        @test @inferred(RM._domain_bounds(:integer,lo,hi)) === (Float64(lo),Float64(hi))
    end
    @test @inferred(RM._domain_bounds(:binary,-2,3.0f0)) === (0.0,1.0)
    @test_throws ArgumentError RM._domain_bounds(:binary,2,3)
    @test_throws MethodError RM._domain_bounds(:continuous,nothing,10)
end

@testset "Helper and functional pairs match DSL literals in build and updates" begin
    for T in (Int,Float32), Backend in (ClpBackend,HiGHSBackend)
        data=(n=4,)
        points=((lo=T(0),hi=T(5),cost=1.0), (lo=T(2),hi=T(4),cost=-1.0),
                (lo=T(0),hi=T(5),cost=1.0))
        sources=(ConstData(HelperBounds(data),data), ConstData(LiteralBounds(data),data),
                 functional(FullBounds(),4), functional(SeparateBounds(),4))
        for source in sources
            lp=LinearProblem(source,first(points); backend=Backend(),safe_mode=true)
            for point in points
                expected=point.cost>0 ? Float64(point.hi) : Float64(point.lo)
                solution=@inferred solve(lp,point)
                @test objective_value(solution)==4 * point.cost * expected
                @test all(k->primal(solution,:x,k)==expected,1:4)
                @test eltype(lp.workspace.lower_values)===Float64
            end
            @test allocated_sweep(lp,points)==0
            clone=similar(lp)
            @test objective_value(clone,points[2])==-8.0
            rebuild!(lp;parameters=points[2])
            @test objective_value(lp,points[2])==-8.0
            @test objective_value(clone,points[1])==20.0
        end
    end
end

@testset "Compact callback return types remain public and numeric storage is normalized" begin
    for T in (Int,Float32), Backend in (ClpBackend,HiGHSBackend)
        point=(cost=1.0,)
        compact=CompactBounds(T(0),T(5))
        @test RM.variable_bounds(compact) isa Tuple{T,T}
        @test RM.variable_bounds(compact,nothing,point,1) isa Tuple{T,T}
        @test RM.variable_bounds(KeyBounds{T}(),1) isa Tuple{T,T}
        for (variable,expected) in ((compact,20.0),(KeyBounds{T}(),26.0))
            lp=LinearProblem(functional(variable,4),point; backend=Backend())
            @test objective_value(lp,point)==expected
            @test objective_value(lp,(cost=2.0,))==2expected
        end
    end
end

@testset "Infinite and mixed pairs preserve Float64 conversion semantics" begin
    data=(n=1,)
    for Backend in (ClpBackend,HiGHSBackend)
        for point in ((lo=-Inf32,hi=3.25f0,cost=1.0), (lo=0.0f0,hi=Inf32,cost=-1.0),
                      (lo=0,hi=Float32(1.2),cost=1.0), (lo=1//3,hi=7//2,cost=1.0))
            lp=LinearProblem(ConstData(HelperBounds(data),data),point; backend=Backend())
            expected=point.cost>0 ? point.cost*Float64(point.hi) : point.cost*Float64(point.lo)
            @test objective_value(lp,point)≈expected
        end
    end
end

@testset "Integer and binary pairs retain domain checks and recovery" begin
    for T in (Int,Float32)
        points=((lo=T(-2),hi=T(3),cost=1.0), (lo=T(-2),hi=T(0),cost=1.0))
        lp=LinearProblem(functional(FullBounds(),1;domain=:binary),points[1];safe_mode=true)
        solution=solve(lp,points[1])
        @test objective_value(solution)==1.0
        @test objective_value(lp,points[2])==0.0
        @test objective_value(lp,points[1])==1.0
        @test_throws ArgumentError solve(lp,(lo=T(2),hi=T(3),cost=1.0))
        @test !lp.built
        @test_throws ErrorException solution.objective_value
        @test objective_value(lp,points[1])==1.0
    end
    p=(lo=1.2f0,hi=4.2f0,cost=-1.0)
    lp=LinearProblem(functional(FullBounds(),1;domain=:integer),p)
    @test objective_value(lp,p)==-2.0
end

@testset "A bound conversion failure invalidates results and recovers" begin
    for Backend in (ClpBackend,HiGHSBackend)
        good=(lo=0,hi=5,cost=1.0,bad=false)
        bad=(lo=0,hi=5,cost=1.0,bad=true)
        lp=LinearProblem(functional(FaultyBounds(),1),good;backend=Backend())
        solution=solve(lp,good)
        @test objective_value(solution)==5.0
        @test_throws MethodError solve(lp,bad)
        @test !lp.built
        @test_throws ErrorException solution.objective_value
        @test objective_value(lp,good)==5.0
        fresh=LinearProblem(functional(FaultyBounds(),1),bad;backend=Backend())
        @test_throws MethodError solve(fresh,bad)
        @test !fresh.built
    end
end

end # module
