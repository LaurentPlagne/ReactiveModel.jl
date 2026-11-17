module ArchitectureTests
using Test
using ReactiveModel

struct MissingRHS <: AbstractConstraint end
ReactiveModel.lhs_terms!(::MissingRHS, emit, data, point, key) = set_coeff!(emit, :x, 1.0)
struct ExplicitZeroRHS <: AbstractConstraint end
ReactiveModel.rhs_value(::ExplicitZeroRHS) = 0.0

@testset "Required RHS and explicit compact RHS methods" begin
    @test_throws MethodError rhs_value(MissingRHS(), nothing, nothing, nothing)
    @test rhs_value(ExplicitZeroRHS(), nothing, nothing, nothing)==0.0
    definition=ProblemDefinition(
        variables=(x=VariableBlock(:x; lower=0.0, upper=10.0),),
        constraints=(floor=ConstraintBlock(:ge, MissingRHS()),),
        objective=Objective(:min, (emit, d, p)->set_coeff!(emit, :x, 1.0)),
    )
    lp=LinearProblem(definition, nothing; safe_mode=true)
    @test_throws MethodError solve(lp, nothing)
    @test !lp.built
end

function native_time_limit(backend::HiGHSBackend)
    value=Ref{Cdouble}()
    status=ccall(
        (:Highs_getDoubleOptionValue, ReactiveModel._HiGHSCAPI.libhighs),
        Cint,
        (Ptr{Cvoid}, Cstring, Ref{Cdouble}),
        backend.inner,
        "time_limit",
        value,
    )
    @assert status==0
    return value[]
end

@model SingleColumn(data, point) begin
    @variable x in [0.0, point.limit]
    @constraint floor x >= point.rhs
    @objective Min x
end

@model ShiftedColumn(data, point) begin
    @variable x in [0.0, point.limit]
    @constraint floor x >= point.rhs + 1.0
    @objective Min 2.0*x
end

client_variables(solution) = solution.variables
client_constraints(solution) = solution.constraints

@noinline function execution_business_terms!(emit, data, point, key)
    point.fail && throw(DomainError(point.rhs, "business callback probe"))
    set_coeff!(emit, :x, data.factor)
    return nothing
end

function execution_sweep!(problem, factor)
    for i in 1:8
        rhs=4.0+i/4
        objective_value(problem, (rhs=rhs, fail=false))==rhs/factor || return false
    end
    return true
end

function counted_clone_rhs(data, point)
    data.evaluations[] += 1
    point.fail && throw(DomainError(point.rhs, "clone assembly probe"))
    return data.factor * point.rhs
end

@model CloneAssembly(data, point) begin
    @variable x in [0.0, 100.0]
    @constraint floor x >= counted_clone_rhs(data, point)
    @objective Min x
end

@testset "Clones with replacement data assemble once without evaluating the source" begin
    p=(rhs=2.0,fail=false)
    source_data=(factor=2.0,evaluations=Ref(0))
    source=LinearProblem(ConstData(CloneAssembly(source_data),source_data),p;backend=ClpBackend())
    @test objective_value(source,p)==4.0
    @test source_data.evaluations[]==1
    replacement=(factor=3.0,evaluations=Ref(0))
    cloned=@inferred similar(source;data=replacement)
    @test source_data.evaluations[]==1
    @test replacement.evaluations[]==1
    @test cloned.const_data.data===replacement
    @test cloned.const_data.definition===source.const_data.definition
    @test cloned.backend.inner!=source.backend.inner
    @test objective_value(cloned,p)==6.0
    @test replacement.evaluations[]==1
    @test objective_value(source,p)==4.0
    @test_throws ArgumentError similar(source;data=(factor=3,evaluations=Ref(0)))
    @test source.numerical_updates.argument.value===nothing
    @test_throws DomainError similar(source,(rhs=2.0,fail=true);data=replacement)
    @test source.numerical_updates.argument.value===nothing
    @test objective_value(source,p)==4.0
    @test source_data.evaluations[]==1
    switched_data=(factor=4.0,evaluations=Ref(0))
    switched=similar(source;data=switched_data,backend=HiGHSBackend())
    @test objective_value(switched,p)==8.0
    @test switched_data.evaluations[]==1
    @test source_data.evaluations[]==1

    unbuilt=LinearProblem(ConstData(CloneAssembly(source_data),source_data),p;backend=ClpBackend())
    deferred_data=(factor=5.0,evaluations=Ref(0))
    deferred=@inferred similar(unbuilt;data=deferred_data)
    @test !deferred.built
    @test deferred_data.evaluations[]==0
    @test objective_value(deferred,p)==10.0
    @test deferred_data.evaluations[]==1
    @test source_data.evaluations[]==1
end

@testset "Execution handles retain callbacks, exceptions and independent clones" begin
    data=(factor=2.0,)
    p=(rhs=4.0, fail=false)
    definition=ProblemDefinition(
        variables=(x=VariableBlock(:x; lower=0.0, upper=10.0),),
        constraints=(floor=ConstraintBlock(
            :floor, (0,), :ge, execution_business_terms!, (d, p, k)->p.rhs,
        ),),
        objective=Objective(:min, (emit, d, p)->set_coeff!(emit, :x, 1.0)),
    )
    source=ConstData(definition, data)
    lp=LinearProblem(source, p; backend=ClpBackend(), safe_mode=true)
    @test lp.const_data===source
    GC.gc()
    solution=@inferred solve(lp, p)
    @test objective_value(solution)==2.0
    @test lp.workspace.result_layout!==nothing
    try
        objective_value(lp, (rhs=5.0, fail=true))
        @test false
    catch error
        @test error isa DomainError
        @test occursin("execution_business_terms!", sprint(Base.show_backtrace, catch_backtrace()))
    end
    @test !lp.built
    @test lp.numerical_updates.argument.value===nothing
    @test_throws ErrorException primal(solution, :x)
    @test objective_value(lp, p)==2.0

    cloned=@inferred similar(lp)
    @test cloned.const_data===lp.const_data
    @test cloned.numerical_updates.payload!==lp.numerical_updates.payload
    @test cloned.numerical_updates.state!==lp.numerical_updates.state
    @test cloned.numerical_updates.argument!==lp.numerical_updates.argument
    @test objective_value(cloned, p)==2.0
    @test objective_value(lp, (rhs=6.0, fail=false))==3.0
    @test objective_value(cloned, p)==2.0

    rebuild!(lp, (factor=4.0,))
    @test lp.const_data.data.factor==4.0
    @test objective_value(lp, p)==1.0
    @test objective_value(cloned, p)==2.0
    current=solve(lp, p)
    # A rejected schema replacement must leave the old model and view usable.
    @test_throws ArgumentError rebuild!(lp, ConstData(SingleColumn(nothing), data))
    @test objective_value(current)==1.0
    @test lp.numerical_updates.argument.value===nothing
    @test objective_value(lp, p)==1.0

    switched=similar(lp; backend=HiGHSBackend())
    @test switched.backend isa HiGHSBackend
    @test objective_value(switched, p)==1.0
    first_task=Threads.@spawn execution_sweep!(lp, 4.0)
    second_task=Threads.@spawn execution_sweep!(cloned, 2.0)
    @test all(fetch, (first_task, second_task))
end

@testset "Result and workspace types are independent of mathematical formulas" begin
    p=(limit=10.0, rhs=2.0)
    original=LinearProblem(ConstData(SingleColumn(nothing)), p; backend=ClpBackend(), safe_mode=true)
    shifted=LinearProblem(ConstData(ShiftedColumn(nothing)), p; backend=ClpBackend(), safe_mode=true)
    @test typeof(original.definition)!==typeof(shifted.definition)
    @test typeof(original)===typeof(shifted)
    @test typeof(original.workspace)===typeof(shifted.workspace)

    first_solution=@inferred solve(original, p)
    second_solution=@inferred solve(shifted, p)
    @test typeof(first_solution)===typeof(second_solution)
    @test objective_value(first_solution)==2.0
    @test objective_value(second_solution)==6.0
    @test @inferred(primal(first_solution, :x))==2.0
    @test @inferred(dual(second_solution, :floor))==2.0
    @test @inferred(activity(second_solution, :floor))==3.0
    @test getfield(first_solution, :variables_cache)===nothing
    @test getfield(second_solution, :constraints_cache)===nothing

    variables=@inferred client_variables(first_solution)
    constraints=@inferred client_constraints(second_solution)
    @test variables.x[0]==2.0
    @test constraints.floor.activities==[3.0]
    @test solve(original, (limit=10.0, rhs=4.0))===first_solution
    @test solve(shifted, (limit=10.0, rhs=4.0))===second_solution
    @test primal(first_solution, :x)==4.0
    @test activity(second_solution, :floor)==5.0
    @test variables.x[0]==2.0
    @test constraints.floor.activities==[3.0]
    @test first_solution.variables.x[0]==4.0
    @test second_solution.constraints.floor.activities==[5.0]

    cloned=similar(shifted)
    @test typeof(cloned.workspace)===typeof(shifted.workspace)
    clone_solution=solve(cloned, p)
    @test typeof(clone_solution)===typeof(second_solution)
    @test objective_value(clone_solution)==6.0
    @test primal(second_solution, :x)==5.0
    rebuild!(shifted)
    @test_throws ErrorException primal(second_solution, :x)
    @test solve(shifted, p)===second_solution
    @test objective_value(second_solution)==6.0
end

@testset "Clones retain configured native options" begin
    p=(limit=10.0, rhs=2.0)
    lp=LinearProblem(
        ConstData(SingleColumn(nothing)),
        p;
        backend=HiGHSBackend(time_limit=7.0, presolve="off"),
        safe_mode=true,
    )
    @test objective_value(lp, p)==2.0
    cloned=similar(lp)
    @test native_time_limit(lp.backend)==7.0
    @test native_time_limit(cloned.backend)==7.0
    @test cloned.backend.inner!=lp.backend.inner
    @test cloned.backend.options!==lp.backend.options
    @test objective_value(cloned, p)==2.0
    ReactiveModel._highs_set_option!(cloned.backend, "time_limit", 8.0)
    @test native_time_limit(lp.backend)==7.0
    @test native_time_limit(cloned.backend)==8.0
end

requires_integer(::Int) = nothing
client_failure(problem) = requires_integer(problem)

@testset "Compact client diagnostics retain opt-in complete types" begin
    problem=LinearProblem(ConstData(SingleColumn(nothing)), (limit=10.0, rhs=2.0))
    compact=sprint(show, typeof(problem);context=:compact=>true)
    complete=ReactiveModel.full_type_string(problem)
    @test occursin("LinearProblem{P=", compact)
    @test !occursin("VariableUpdater", complete)
    @test !occursin("ProblemDefinition", complete)
    @test problem.numerical_updates.model_name===:SingleColumn
    @test length(complete)>length(compact)
    @test occursin(compact, sprint(show, Tuple{typeof(problem),Int};context=:compact=>true))
    @test !occursin("VariableUpdater", repr(problem))
    # The trim verifier prints free TypeVars inside bounds of inferred types.
    # Those values do not necessarily bind a `::Type{T}` static parameter.
    free_workspace = Workspace{TypeVar(:P)}
    @test sprint(show, free_workspace; context=:compact=>true) == sprint(show, free_workspace)
    @test ReactiveModel.full_type_string(free_workspace) == sprint(show, free_workspace)
    try
        client_failure(problem)
        @test false
    catch error
        @test error isa MethodError
        short_error=sprint(showerror, error)
        io=IOBuffer()
        showerror(IOContext(io, :reactivemodel_full_types=>true), error)
        full_error=String(take!(io))
        @test occursin("LinearProblem{P=", short_error)
        @test length(full_error)>length(short_error)
        trace=sprint(Base.show_backtrace, catch_backtrace())
        @test occursin("LinearProblem{P=", trace)
    end
end

# Result storage can differ from the included native backends. No forwarded
# properties or conventional result fields are provided by this adapter.
struct ResultAdapter{B<:AbstractBackend} <: AbstractBackend
    solver::B
end
ReactiveModel.backend_status(b::ResultAdapter) = ReactiveModel.backend_status(b.solver)
function ReactiveModel.backend_objective_value(b::ResultAdapter)
    return ReactiveModel.backend_objective_value(b.solver)
end
function ReactiveModel.backend_primal_value(b::ResultAdapter, i::Int)
    return ReactiveModel.backend_primal_value(b.solver, i)
end
function ReactiveModel.backend_dual_value(b::ResultAdapter, i::Int)
    return ReactiveModel.backend_dual_value(b.solver, i)
end
for operation in (
    :initialize!,
    :backend_set_verbose_flag!,
    :backend_solve!,
    :backend_change_col_bounds!,
    :update_variable_block_bounds!,
    :update_constraint_rhs!,
    :update_constraint_variable_coefficients!,
    :backend_update_objective_coefficients!,
)
    @eval ReactiveModel.$operation(b::ResultAdapter, args...; kwargs...) =
        ReactiveModel.$operation(b.solver, args...; kwargs...)
end
ReactiveModel.backend_has_duals(b::ResultAdapter) = backend_has_duals(b.solver)
function ReactiveModel.update_constraint_block_rhs!(
    b::ResultAdapter, block::ReactiveModel._BackendConstraintBlock, offset::Int
)
    return ReactiveModel.update_constraint_block_rhs!(b.solver, block, offset)
end

@testset "Dispatched result protocol supports independent storage" begin
    p=(limit=10.0, rhs=2.0)
    lp=LinearProblem(
        ConstData(SingleColumn(nothing)), p; backend=ResultAdapter(ClpBackend()), safe_mode=true
    )
    solution=solve(lp, p)
    @test objective_value(solution)==2.0
    @test primal(solution, :x)==2.0
    @test dual(solution, :floor)==1.0
    @test activity(solution, :floor)==2.0
    @test solution.variables.x[0]==2.0
    @test solution.constraints.floor[0]==1.0
    @test objective_value(lp, (limit=9.0, rhs=3.0))==3.0
    @test primal(solution, :x)==3.0
end

function many_families(n)
    names=ntuple(i->Symbol(:x, i), n)
    rownames=ntuple(i->Symbol(:floor, i), n)
    variables=NamedTuple{names}(
        map(names) do name
            return VariableBlock(name; lower=0.0, upper=(d, p, k)->p.limit)
        end,
    )
    constraints=NamedTuple{rownames}(
        map(names, rownames) do name, row
            return ConstraintBlock(
                row, (nothing,), :ge, (emit, d, p, k)->set_coeff!(emit, name, 1.0), (d, p, k)->p.rhs
            )
        end,
    )
    objective=Objective(:min, (emit, d, p)->foreach(name->set_coeff!(emit, name, 1.0), d.names))
    data=(names=collect(names),)
    return ConstData(ProblemDefinition(; variables, constraints, objective), data)
end
function sweep!(lp)
    for rhs in (1.0, 2.0, 3.0)
        objective_value(lp, (limit=10.0, rhs=rhs))
    end
    return nothing
end
plan_widths(::ReactiveModel._EmptyFamilyPlan) = Int[]
plan_widths(plan::ReactiveModel._FamilyPlan) = [length(plan.head); plan_widths(plan.tail)]

@testset "Bounded plans and key-type indices retain zero-allocation updates" begin
    small=LinearProblem(many_families(2), (limit=10.0, rhs=1.0); backend=ClpBackend())
    larger=LinearProblem(many_families(17), (limit=10.0, rhs=1.0); backend=ClpBackend())
    @test typeof(small)===typeof(larger)
    @test typeof(small.workspace)===typeof(larger.workspace)
    @test typeof(solve(small, (limit=10.0, rhs=1.0)))===typeof(solve(larger, (limit=10.0, rhs=1.0)))
    @test typeof(small.variable_indices)===typeof(larger.variable_indices)
    @test typeof(small.update_scratch)===typeof(larger.update_scratch)
    @test length(larger.variable_indices.groups)==1
    @test length(larger.variable_indices)==17
    implementation=larger.numerical_updates.payload[].updates
    @test plan_widths(implementation.variables)==[8, 8, 1]
    @test plan_widths(implementation.constraints)==[8, 8, 1]
    sweep!(larger)
    @test objective_value(larger, (limit=10.0, rhs=2.0))==34.0
    sweep!(larger)
    @test @allocated(sweep!(larger))==0
    @test primal(solve(larger, (limit=10.0, rhs=2.0)), :x17)==2.0
end
end
