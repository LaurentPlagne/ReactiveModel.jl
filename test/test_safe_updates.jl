module SafeUpdateTests
using Test
using ReactiveModel

# Business functions receive the study and complete point, including branches.
capacity(d,p) = d.capacity * clamp(p.availability,0.0,1.0) *
                (p.temperature > 25.0 ? 1.0-d.derating*(p.temperature-25.0) : 1.0)
conversion(d,p) = d.efficiency * (p.boost ? 1.1 : 1.0) / (1.0+0.001*(p.temperature-20.0)^2)
demand(d,p) = d.load * p.load + max(0.0,p.temperature-25.0)*d.cooling
cost(d,p) = d.fuel*p.price
offset(d,p) = p.boost ? d.startup : 0.0

@model WrongConditions(d,p) begin
    @variable x in [0.0, capacity(d,p)] update_when false
    @variable backup in [0.0, 100.0] update_when false
    @constraint balance conversion(d,p)*x + (p.boost ? 1.0 : 0.0)*backup >= demand(d,p) update_when false
    @objective Min cost(d,p)*x + 100.0*backup + offset(d,p) update_when false
end

@model BatchedConditions(d,p) begin
    @variable x[i in 1:5] in [0.0,p.limit] update_when false
    @constraint floor[i in 1:5] p.a*x[i] >= p.rhs+i update_when false
    @objective Min sum(p.cost*x[i] for i in 1:5) + p.offset update_when false
end

@model RedundantConditions(d,p) begin
    @variable x in [0.0,d.capacity] update_when true
    @constraint floor 1.0*x >= 2.0 update_when true
    @objective Min x + d.startup update_when true
end

findcheck(report,family,component) = only(c for c in report if c.family===family && c.component===component)

@testset "Safe updates correct every skipped numerical component" begin
    data=(capacity=30.0,derating=0.01,efficiency=0.9,load=10.0,cooling=0.1,fuel=2.0,startup=3.0)
    points=((availability=1.0,temperature=20.0,boost=false,load=0.7,price=1.0),
            (availability=0.9,temperature=30.0,boost=true,load=1.1,price=1.3),
            (availability=0.8,temperature=25.0,boost=false,load=0.9,price=0.9))
    for backend in (HiGHSBackend,ClpBackend)
        study=ConstData(WrongConditions(data),data)
        checked=LinearProblem(study,points[1];backend=backend(),safe_mode=true)
        @test isempty(update_report(checked))
        for p in points
            value=objective_value(checked,p)
            fresh=LinearProblem(study,p;backend=backend())
            @test value ≈ objective_value(fresh,p)
            @test checked.upper_values[1] ≈ capacity(data,p)
            @test checked.row_coefficients[1][1] ≈ conversion(data,p)
            @test get(checked.row_coefficients[1],2,0.0) == (p.boost ? 1.0 : 0.0)
            @test checked.row_rhs[1] ≈ demand(data,p)
            @test checked.objective_coefficients[1] ≈ cost(data,p)
            @test checked.objective_offset[] ≈ offset(data,p)
        end
        report=update_report(checked)
        for (family,component) in ((:x,:bounds),(:balance,:coefficients),(:balance,:rhs),(:objective,:objective))
            c=findcheck(report,family,component)
            @test c.checks==2
            @test c.missed==2
            @test c.changes==2
        end
        @test findcheck(report,:backup,:bounds).changes==0
        # Reports are snapshots; cloning inherits the mode but owns its report.
        empty!(report)
        @test !isempty(update_report(checked))
        cloned=similar(checked)
        @test cloned.update_safety.enabled
        @test isempty(update_report(cloned))
        @test objective_value(cloned,points[1]) ≈ objective_value(LinearProblem(study,points[1];backend=backend()),points[1])
        @test length(update_report(checked))==5
        fast=similar(checked;safe_mode=false)
        @test !fast.update_safety.enabled
        @test isempty(update_report(fast))
        @test reset_update_report!(checked) === checked
        @test isempty(update_report(checked))
        # The same-point solve shortcut must not bypass checked evaluation.
        sol=solve(checked,points[end])
        @test all(c.checks==1 for c in update_report(checked))
        @test solve(checked,points[end]) === sol
        @test all(c.checks==2 for c in update_report(checked))
        # A numerical rebuild preserves the mode and existing observations.
        count=length(update_report(checked))
        rebuild!(checked)
        @test checked.update_safety.enabled
        @test length(update_report(checked))==count
        @test objective_value(checked,points[1]) ≈ objective_value(LinearProblem(study,points[1];backend=backend()),points[1])
    end
end

@testset "Safe batched changes, offsets, infeasibility and recovery" begin
    for backend in (HiGHSBackend,ClpBackend)
        p=(limit=20.0,a=1.0,rhs=1.0,cost=2.0,offset=4.0)
        study=ConstData(BatchedConditions(nothing))
        lp=LinearProblem(study,p;backend=backend(),safe_mode=true)
        objective_value(lp,p)
        for q in ((limit=15.0,a=2.0,rhs=2.0,cost=3.0,offset=5.0),
                  (limit=15.0,a=2.0,rhs=2.0,cost=3.0,offset=8.0))
            @test objective_value(lp,q) ≈ objective_value(LinearProblem(study,q;backend=backend()),q)
            @test all(==(q.limit),lp.upper_values)
            @test lp.row_rhs ≈ [q.rhs+i for i in 1:5]
        end
        @test findcheck(update_report(lp),:objective,:objective).missed==2
        sol=solve(lp,lp.parameter_snapshot)
        bad=(limit=1.0,a=1.0,rhs=10.0,cost=2.0,offset=0.0)
        @test_throws SolveError solve(lp,bad)
        @test_throws Exception primal(sol,:x,1)
        @test_throws SolveError solve(lp,bad)
        @test objective_value(lp,p) ≈ 44.0
    end
end

@testset "Unnecessary means observed numerical invariance" begin
    data=(capacity=10.0,startup=3.0)
    lp=LinearProblem(ConstData(RedundantConditions(data),data),(step=1,);safe_mode=true)
    @test objective_value(lp,(step=1,))==5.0
    @test objective_value(lp,(step=2,))==5.0
    report=update_report(lp)
    @test length(report)==4
    @test all(c.checks==1 && c.changes==0 && c.missed==0 && c.unnecessary==1 for c in report)
    fast=LinearProblem(lp.const_data,(step=1,))
    @test objective_value(fast,(step=1,))==5.0
    @test objective_value(fast,(step=2,))==5.0
    @test isempty(update_report(fast))
end

@testset "Functional callbacks and recovery from an evaluation exception" begin
    definition=ProblemDefinition(
        variables=(x=VariableBlock(:x;lower=(d,p,k)->p.lower,
            upper=(d,p,k)->p.limit,has_changed=(old,new)->false),),
        constraints=(floor=ConstraintBlock(:floor,(nothing,),:ge,
            (emit,d,p,k)->set_coeff!(emit,:x,p.a),
            (d,p,k)->p.rhs < 0 ? throw(ArgumentError("negative demand")) : p.rhs;
            has_changed=(old,new)->false),),
        objective=Objective(:min,(emit,d,p)->set_coeff!(emit,:x,p.cost);
            has_changed=(old,new)->false))
    p=(lower=0.0,limit=10.0,a=1.0,rhs=2.0,cost=3.0)
    q=(lower=4.0,limit=12.0,a=2.0,rhs=6.0,cost=5.0)
    for backend in (HiGHSBackend,ClpBackend)
        lp=LinearProblem(definition,p;backend=backend(),safe_mode=true)
        sol=@inferred solve(lp,p)
        @test objective_value(sol)==6.0
        @test objective_value(lp,q)==20.0
        @test lp.lower_values==[4.0]
        @test all(c.missed==1 for c in @inferred(update_report(lp)))
        @test_throws ArgumentError solve(lp,merge(q,(rhs=-1.0,)))
        @test !lp.built
        @test_throws ErrorException objective_value(sol)
        @test objective_value(lp,p)==6.0
    end
end
end
