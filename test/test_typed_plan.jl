using Test
using ReactiveModel
include("trim/workload.jl")
using .ReactiveModelTrimWorkload

collection_scalars(solution) = solution.variables.label[:first] + solution.variables.x[1] + solution.constraints.floor[1]
@noinline function retain_one_result_family(problem)
    solution=solve(problem, (a=1.0, rhs=1.0, cost=2.0, extra=0.0))
    collection=solution.variables
    unrelated=WeakRef(collection.label.values)
    retained=collection.x
    objective_value(problem, (a=1.0, rhs=2.0, cost=2.0, extra=0.0))
    return retained, unrelated
end

@testset "Concrete update plan and result types" begin
    small = @inferred ReactiveModelTrimWorkload.lp_problem(1)
    larger = @inferred ReactiveModelTrimWorkload.lp_problem(5)
    # Runtime dimensions do not create a new compiled model family.
    @test typeof(small) === typeof(larger)
    @test isconcretetype(typeof(small.numerical_updates))
    @test isconcretetype(fieldtype(typeof(small.variable_indices),:groups))
    implementation = small.numerical_updates.payload[].updates
    @test isconcretetype(fieldtype(typeof(implementation), :variables))
    @test implementation.variables isa ReactiveModel._AbstractFamilyPlan
    @test implementation.constraints isa ReactiveModel._AbstractFamilyPlan
    @test_throws ArgumentError merge(small.definition.variables, small.definition.variables)
    @test_throws ArgumentError merge(small.definition.constraints, small.definition.constraints)
    clone = similar(small)
    @test typeof(clone) === typeof(small)
    @test objective_value(clone, (a=1.0, rhs=1.0, cost=2.0, extra=0.0)) == 2.0

    sol = @inferred solve(larger, (a=1.0, rhs=1.0, cost=2.0, extra=0.0))
    @test @inferred(primal(sol, :label, :first)) == 0.0
    @test @inferred(dual(sol, :floor, 1)) == 2.0
    @test @inferred(activity(sol, :floor, 1)) == 1.0
    @test getfield(sol, :variables_cache) === nothing
    @test getfield(sol, :constraints_cache) === nothing
    @test @inferred(sol.variables.label[:first]) == 0.0
    @test eltype(sol.variables.label.keys) === Symbol
    @test eltype(sol.variables.pair.keys) === Tuple{Symbol,Int}
    @test eltype(sol.variables.x.keys) === Int
    @test eltype(sol.constraints.named.keys) === Symbol
    @test keys(sol.variables) == [:x, :anchor, :label, :pair, :hole]
    @test [block.name for block in sol.variables.blocks] == keys(sol.variables)
    @test sol.variables.by_name[:label][:first] == 0.0
    @test :keys in propertynames(sol.variables.label)
    @test !occursin("SolutionCollection", repr(sol.variables.label))
    @test @inferred(collection_scalars(sol)) == 3.0
    collection_scalars(sol)
    @test @allocated(collection_scalars(sol)) == 0
    # Scalar keys normalize consistently before and after collection access.
    before = dual(sol, :zero)
    sol.constraints
    @test dual(sol, :zero) == before
    retained, unrelated=retain_one_result_family(ReactiveModelTrimWorkload.lp_problem(2))
    GC.gc()
    @test retained[1] == 1.0
    @test unrelated.value === nothing

    @test ReactiveModelTrimWorkload.verify(1)
    @test ReactiveModelTrimWorkload.verify(5)
    @test ReactiveModelTrimWorkload.verify_math(9)
    @test ReactiveModelTrimWorkload.verify_lp(5; backend=ClpBackend())
end

module GroupedPlanConstructionValidation
using Test
using ReactiveModel

struct GroupVariable{N} <: AbstractVariable end
struct GroupRow{N} <: AbstractConstraint end
ReactiveModel.lower_bound(::GroupVariable{N}, data, point, key) where {N} = 0.0
ReactiveModel.upper_bound(::GroupVariable{N}, data, point, key) where {N} = 100.0
function ReactiveModel.lhs_terms!(::GroupRow{N}, emit, data, point, key) where {N}
    return set_coeff!(emit, Symbol(:x, N), key, point.a)
end
ReactiveModel.rhs_value(::GroupRow{N}, data, point, key) where {N} =
    data.scale[1] * point.rhs

function grouped_problem(nfamilies=17; scale=1.0)
    names = ntuple(i -> Symbol(:x, i), nfamilies)
    rownames = ntuple(i -> Symbol(:floor, i), nfamilies)
    # Empty families straddle group boundaries while family 17 stays nonempty.
    sizes = [1,2,3,1,2,3,1,0,0,2,3,1,2,3,1,0,3][1:nfamilies]
    variables = NamedTuple{names}(ntuple(nfamilies) do i
        VariableBlock(names[i], 1:sizes[i], GroupVariable{i}())
    end)
    constraints = NamedTuple{rownames}(ntuple(nfamilies) do i
        ConstraintBlock(rownames[i], 1:sizes[i], :ge, GroupRow{i}())
    end)
    objective = Objective(:min, (emit,d,p) -> begin
        for (i,name) in enumerate(d.names), k in 1:d.sizes[i]
            set_coeff!(emit, name, k, 1.0)
        end
    end)
    data = (names=Symbol[names...], sizes=sizes, scale=[scale])
    definition = ProblemDefinition(; variables, constraints, objective)
    return LinearProblem(ConstData(definition,data), (a=1.0,rhs=1.0);
        backend=HiGHSBackend(output_flag=false)), data
end

widths(::ReactiveModel._EmptyFamilyPlan) = Int[]
widths(plan::ReactiveModel._FamilyPlan) = [length(plan.head); widths(plan.tail)]

function check_offsets(lp, sizes)
    plan = lp.numerical_updates.payload[].updates
    nfamilies = length(sizes)
    expected_widths = fill(8, nfamilies ÷ 8)
    nfamilies % 8 == 0 || push!(expected_widths, nfamilies % 8)
    @test widths(plan.variables) == widths(plan.constraints) == expected_widths
    indices, ranges = Int[], UnitRange{Int}[]
    ReactiveModel._foreach_family(plan.variables) do updater
        push!(indices, updater.block_index)
        push!(ranges, updater.columns)
    end
    starts = [1 + sum(sizes[1:i-1]) for i in eachindex(sizes)]
    @test indices == collect(1:nfamilies)
    @test ranges == [start:(start+n-1) for (start,n) in zip(starts,sizes)]
    empty!(indices)
    rows = Int[]
    ReactiveModel._foreach_family(plan.constraints) do updater
        push!(indices, updater.block_index)
        push!(rows, updater.first_row)
    end
    @test indices == collect(1:nfamilies)
    @test rows == starts
end

@testset "Grouped construction preserves family indices and empty offsets" begin
    for nfamilies in (0,8,9,16,17)
        boundary, boundary_data = grouped_problem(nfamilies)
        @test !boundary.built
        @test length(Set(typeof(block.variable) for block in boundary.definition.variables.blocks)) == nfamilies
        @test length(Set(typeof(block.constraint) for block in boundary.definition.constraints.blocks)) == nfamilies
        check_offsets(boundary, boundary_data.sizes)
    end
    lp, data = grouped_problem()
    p, q = (a=1.0,rhs=1.0), (a=2.0,rhs=3.0)
    count = sum(data.sizes)
    plan_type = typeof(lp.numerical_updates.payload[].updates)
    check_offsets(lp, data.sizes)
    @test objective_value(lp,p) ≈ count
    @test objective_value(lp,q) ≈ 1.5count
    fresh, _ = grouped_problem()
    @test objective_value(lp,q) ≈ objective_value(fresh,q)
    @test primal(solve(lp,q),:x17,3) ≈ 1.5
    check_offsets(lp, data.sizes)
    clone = similar(lp;data=(names=copy(data.names),sizes=copy(data.sizes),scale=[1.0]))
    data.scale[1] = 2.0
    rebuild!(lp)
    @test typeof(lp.numerical_updates.payload[].updates) === plan_type
    check_offsets(lp, data.sizes)
    @test objective_value(lp,q) ≈ 3.0count
    rebuilt_fresh, _ = grouped_problem(;scale=2.0)
    @test objective_value(lp,q) ≈ objective_value(rebuilt_fresh,q)
    @test objective_value(clone,q) ≈ 1.5count
    @test primal(solve(lp,q),:x17,3) ≈ 3.0
end
end
