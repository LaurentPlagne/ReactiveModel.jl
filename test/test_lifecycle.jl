module LifecycleTests
using Test
using ReactiveModel

@noinline function checked_value(data, stage::Symbol, value)
    data.calls[] += 1
    data.failure[] === stage && throw(DomainError(stage, "lifecycle callback probe"))
    return value
end
@noinline function checked_rhs(data, point)
    data.failure[] === :solve_error && throw(SolveError(:other_error, "callback SolveError probe"))
    return checked_value(data, :rhs, data.minimum[] * point.scale)
end
probe_data(value=2.0) = (minimum=Ref(value), failure=Ref(:none), calls=Ref(0))

@model LifecycleModel(data, point) begin
    @variable x in [checked_value(data, :bounds, 0.0), 20.0]
    @constraint floor checked_value(data, :coefficients, 1.0) * x >= checked_rhs(data, point)
    @objective Min checked_value(data, :objective, 1.0) * x
end

function problem(data, point, backend; safe_mode=false)
    return LinearProblem(ConstData(LifecycleModel(data), data), point; backend, safe_mode)
end

live_objective(solution) = solution.objective_value
live_status(solution) = solution.termination_status
scalar_properties_match(solution, value) = live_objective(solution) == value && live_status(solution) === :optimal

function check_stale(solution)
    @test_throws ErrorException solution.objective_value
    @test_throws ErrorException solution.termination_status
    @test_throws ErrorException solution.variables
    @test_throws ErrorException solution.constraints
    @test_throws ErrorException objective_value(solution)
    @test_throws ErrorException termination_status(solution)
    @test_throws ErrorException is_optimal(solution)
    @test_throws ErrorException primal(solution, :x)
    @test_throws ErrorException dual(solution, :floor)
    @test_throws ErrorException activity(solution, :floor)
end

function check_rebuild_required(lp)
    @test !lp.built
    @test lp.parameter_snapshot === nothing
    @test lp.last_status === :not_solved
    @test lp.numerical_updates.argument.value === nothing
end

@testset "Interrupted rebuild invalidates all reads and recovers from callbacks" begin
    p, q = (scale=1.0,), (scale=2.0,)
    for factory in (ClpBackend, HiGHSBackend), stage in (:bounds, :coefficients, :rhs, :objective), next in (p, q)
        data = probe_data()
        lp = problem(data, p, factory())
        solution = solve(lp, p)
        saved_variables, saved_constraints = solution.variables, solution.constraints
        clone = similar(lp; data=probe_data())
        @test objective_value(solution) == 2.0
        @test objective_value(clone, p) == 2.0

        data.minimum[] = 3.0
        data.failure[] = stage
        @test_throws DomainError rebuild!(lp)
        check_rebuild_required(lp)
        check_stale(solution)
        @test saved_variables.x[0] == 2.0
        @test saved_constraints.floor[0] == 1.0
        calls = data.calls[]
        @test_throws DomainError objective_value(lp, p)
        @test data.calls[] > calls
        check_rebuild_required(lp)
        @test objective_value(clone, p) == 2.0

        data.failure[] = :none
        fresh = problem(probe_data(3.0), next, factory())
        @test objective_value(lp, next) == objective_value(fresh, next) == 3.0 * next.scale
        @test solve(lp, next) === solution
        @test solution.objective_value == objective_value(solution)
        @test solution.termination_status === termination_status(solution) === :optimal
        @test primal(solution, :x) == 3.0 * next.scale
        @test objective_value(clone, p) == 2.0
    end
end

# Mutate the native model and then fail too: clearing only RM's buffers is not
# enough to prove recovery after an initialization error.
mutable struct InitializationProbe{B<:AbstractBackend} <: AbstractBackend
    solver::B
    failure::Symbol
end
function ReactiveModel.initialize!(backend::InitializationProbe, args...)
    backend.failure === :before && throw(ErrorException("before initialize probe"))
    ReactiveModel.initialize!(backend.solver, args...)
    backend.failure === :after && throw(ErrorException("after initialize probe"))
    return nothing
end
for operation in (:backend_set_verbose_flag!, :backend_solve!, :backend_change_col_bounds!,
    :backend_change_cols_bounds!, :update_variable_block_bounds!, :update_constraint_rhs!,
    :update_constraint_block_rhs!, :update_constraint_variable_coefficients!,
    :backend_update_objective_coefficients!, :backend_status, :backend_objective_value,
    :backend_primal_value, :backend_dual_value, :backend_has_duals)
    @eval ReactiveModel.$operation(b::InitializationProbe, args...; kwargs...) =
        ReactiveModel.$operation(b.solver, args...; kwargs...)
end

@testset "Interrupted native initialization requires a fresh assembly" begin
    p, q = (scale=1.0,), (scale=2.0,)
    for factory in (ClpBackend, HiGHSBackend), stage in (:before, :after), next in (p, q)
        data = probe_data()
        backend = InitializationProbe(factory(), :none)
        lp = problem(data, p, backend)
        solution = solve(lp, p)
        @test solution.objective_value == 2.0
        data.minimum[] = 3.0
        backend.failure = stage
        @test_throws ErrorException rebuild!(lp)
        check_rebuild_required(lp)
        check_stale(solution)
        @test_throws ErrorException solve(lp, p)
        check_rebuild_required(lp)
        backend.failure = :none
        @test objective_value(lp, next) == objective_value(problem(probe_data(3.0), next, factory()), next)
        @test solve(lp, next) === solution
        @test solution.objective_value == 3.0 * next.scale
    end
end

@testset "Public scalar properties follow failure and successful live refresh" begin
    p, bad = (scale=1.0,), (scale=100.0,)
    for factory in (ClpBackend, HiGHSBackend), safe_mode in (false, true)
        lp = problem(probe_data(), p, factory(); safe_mode)
        solution = solve(lp, p)
        @test (@inferred live_objective(solution)) == (@inferred objective_value(solution)) == 2.0
        @test (@inferred live_status(solution)) === (@inferred termination_status(solution)) === :optimal
        @test scalar_properties_match(solution, 2.0)
        @test @allocated(scalar_properties_match(solution, 2.0)) == 0
        @test_throws SolveError solve(lp, bad)
        @test lp.built
        @test lp.parameter_snapshot == Parameters(bad)
        @test lp.last_status === :infeasible
        check_stale(solution)
        @test_throws SolveError objective_value(lp, bad)
        @test lp.built
        @test solve(lp, p) === solution
        @test solution.objective_value == 2.0
        @test solution.termination_status === :optimal
    end
end

@testset "A callback SolveError has the same recovery policy as other exceptions" begin
    p, q = (scale=1.0,), (scale=2.0,)
    for factory in (ClpBackend, HiGHSBackend), entry in (:rebuild, :update)
        data = probe_data()
        lp = problem(data, p, factory())
        solution = solve(lp, p)
        data.failure[] = :solve_error
        if entry === :rebuild
            @test_throws SolveError rebuild!(lp; parameters=q)
        else
            @test_throws SolveError objective_value(lp, q)
        end
        check_rebuild_required(lp)
        check_stale(solution)
        data.failure[] = :none
        @test objective_value(lp, q) == 4.0
        @test solve(lp, q) === solution
    end
end

@testset "Deferred rebuild and rejected schema preserve their distinct contracts" begin
    p = (scale=1.0,)
    for factory in (ClpBackend, HiGHSBackend)
        data = probe_data()
        lp = problem(data, p, factory())
        solution = solve(lp, p)
        data.minimum[] = 3.0
        rebuild!(lp; parameters=nothing)
        check_rebuild_required(lp)
        check_stale(solution)
        clone = similar(lp)
        @test !clone.built
        @test objective_value(clone, p) == 3.0
        @test solve(lp, p) === solution
        @test solution.objective_value == 3.0
        old_rows = copy(lp.row_rhs)
        other_definition = ProblemDefinition(
            variables=(y=VariableBlock(:y; lower=0.0, upper=10.0),),
            constraints=NamedTuple(),
            objective=Objective(:min, (emit, d, point)->set_coeff!(emit, :y, 1.0)),
        )
        @test_throws ArgumentError rebuild!(lp, ConstData(other_definition, data))
        @test lp.built
        @test lp.row_rhs == old_rows
        @test lp.last_status === :optimal
        @test solve(lp, p) === solution
        @test solution.objective_value == 3.0
        @test solution.termination_status === :optimal
        @test lp.numerical_updates.argument.value === nothing
    end
end
end # module
