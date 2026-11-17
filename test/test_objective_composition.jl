module ObjectiveCompositionTests
using Test
using ReactiveModel

@model SharedCost(data, point) begin
    @external_variables x
    @objective Min data.scale[1] * point.a * x + point.offset update_when (:a, :offset)
end

@model SkippedCost(data, point) begin
    @external_variables x
    @objective Min point.a * x + point.offset update_when false
end

# Instrument a concrete child without changing its mathematical contribution.
struct TracedCost{O} <: AbstractObjective
    child::O
    label::Int
end
function ReactiveModel.objective_terms!(cost::TracedCost, emit, data, point)
    push!(data.calls, cost.label)
    objective_terms!(cost.child, emit, data, point)
end
ReactiveModel.has_changed(cost::TracedCost, old, new) = has_changed(cost.child, old, new)

point(; a=2.0, b=3.0, offset=4.0, active=true, fail=false, tick=0) =
    (; a, b, offset, active, fail, tick)
study() = (scale=[2.0], calls=Int[])
function definition(data; skipped=false)
    dsl = skipped ? SkippedCost(data).objective : SharedCost(data).objective
    first = Objective(:min, TracedCost(dsl.objective, 1))
    second = Objective(:min, (emit, d, p) -> begin
        push!(d.calls, 2)
        p.fail && throw(DomainError(p.b, "second objective fragment"))
        p.active && set_coeff!(emit, :x, p.b)
        ReactiveModel._add_constant!(emit, 1.0)
    end; has_changed=skipped ? ((old,new)->false) :
        ((old,new)->old.b != new.b || old.active != new.active || old.fail != new.fail))
    return ProblemDefinition(
        variables=(x=VariableBlock(:x; lower=1.0, upper=4.0),),
        constraints=NamedTuple(), objective=compose_objectives(first, second))
end

@testset "Objective composition contract and bounded traversal" begin
    a = Objective(:min, (emit,d,p)->set_coeff!(emit,:x,1.0))
    b = Objective(:max, (emit,d,p)->nothing)
    @test compose_objectives(a) === a
    @test_throws MethodError compose_objectives()
    @test_throws ArgumentError compose_objectives(a,b)
    @test_throws ArgumentError compose_objectives(a,a,b)
    # Nesting and more than eight children exercise both traversal groups.
    nested = compose_objectives(a, compose_objectives(ntuple(_->a, 10)...))
    def = ProblemDefinition(variables=(x=VariableBlock(:x;lower=1.0,upper=4.0),),
        constraints=NamedTuple(),objective=nested)
    lp = LinearProblem(def,(tick=0,))
    @test objective_value(lp,(tick=0,)) == 11.0
    @test objective_value(lp,(tick=1,)) == 11.0
    # Only the last child changes, beyond the first group; early true predicates
    # must short-circuit before the deliberately throwing later predicate.
    never = Objective(:min,(emit,d,p)->nothing;has_changed=(old,new)->false)
    last = Objective(:min,(emit,d,p)->nothing;has_changed=(old,new)->old != new)
    grouped = compose_objectives(ntuple(_->never,9)...,last)
    @test !has_changed(grouped.objective,(tick=0,),(tick=0,))
    @test has_changed(grouped.objective,(tick=0,),(tick=1,))
    bomb = Objective(:min,(emit,d,p)->nothing;
        has_changed=(old,new)->error("predicate must be skipped"))
    @test has_changed(compose_objectives(last,bomb).objective,(tick=0,),(tick=1,))
end

@testset "Complete ordered reemission, cancellations and recovery" begin
    for backend in (HiGHSBackend,ClpBackend)
        data=study()
        p=point()
        def=definition(data)
        lp=LinearProblem(ConstData(def,data),p;backend=backend())
        @test objective_value(lp,p)==12.0 # (2*2 + 3)*1 + 4 + 1
        @test data.calls==[1,2]
        empty!(data.calls)
        q=point(a=3.0)
        @test objective_value(lp,q)==14.0
        @test data.calls==[1,2] # unchanged second child is emitted too
        empty!(data.calls)
        @test objective_value(lp,point(a=3.0,tick=1))==14.0
        @test isempty(data.calls)
        # Delete the second contribution, cancel the entire column, then return.
        for r in (point(active=false),point(b=-4.0,offset=7.0),
                  point(a=0.0,b=0.0,offset=0.0),p)
            freshdata=study()
            fresh=LinearProblem(ConstData(definition(freshdata),freshdata),r;backend=backend())
            @test objective_value(lp,r) ≈ objective_value(fresh,r)
            @test get(lp.objective_coefficients,1,0.0) ==
                data.scale[1]*r.a + (r.active ? r.b : 0.0)
            @test lp.objective_offset[]==r.offset+1.0
        end
        @test_throws DomainError objective_value(lp,point(fail=true))
        @test !lp.built
        empty!(data.calls)
        @test objective_value(lp,p)==12.0
        @test data.calls==[1,2]
        independent=similar(lp;data=study())
        @test independent.const_data.data.scale !== data.scale
        @test independent.const_data.data.calls !== data.calls
        data.scale[1]=4.0
        rebuild!(lp)
        @test objective_value(lp,p)==16.0
        @test objective_value(independent,p)==12.0
        replacement=(scale=[3.0],calls=Int[])
        rebuild!(lp,replacement)
        @test objective_value(lp,p)==14.0
    end
end

@testset "Safe diagnostics describe one aggregate objective" begin
    for backend in (HiGHSBackend,ClpBackend)
        data=study()
        p=point()
        lp=LinearProblem(ConstData(definition(data;skipped=true),data),p;
            backend=backend(),safe_mode=true)
        @test objective_value(lp,p)==10.0
        reset_update_report!(lp)
        empty!(data.calls)
        q=point(a=4.0,b=5.0,offset=8.0)
        @test objective_value(lp,q)==18.0
        @test data.calls==[1,2]
        checks=filter(c->c.component===:objective,update_report(lp))
        @test length(checks)==1
        check=only(checks)
        @test check.family===:objective
        @test check.checks==1 && check.missed==1 && check.changes==1
    end
end
end
