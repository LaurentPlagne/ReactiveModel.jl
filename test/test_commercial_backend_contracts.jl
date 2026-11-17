module CommercialBackendContractTests

using Test
using ReactiveModel
const RM = ReactiveModel

@testset "CPLEX unsupported operations fail without native access" begin
    # Null handles deliberately cannot call a solver: these capability guards
    # are part of RM, independently of extension availability or a license.
    backend = CPLEXBackend(C_NULL, C_NULL, Float64[], Float64[], Float64[],
                           Float64[], false, NaN, :not_solved, Dict{Symbol,Any}())
    variable = RM._BackendVariableBlock(:continuous, 0, [0.0], [1.0], [0.0])
    constraint = RM._BackendConstraint([1], [2.0], 3.0, :le)
    block = RM._BackendConstraintBlock(0, [constraint], [0.0])
    objective = RM._BackendObjectiveCoefficients([1], [2.0])
    operations = (
        () -> backend_change_col_bounds!(backend, 1, 0.0, 2.0),
        () -> backend_change_cols_bounds!(backend, [1], [0.0], [2.0]),
        () -> update_variable_block_bounds!(backend, variable, 0),
        () -> update_constraint_rhs!(backend, constraint, 1),
        () -> update_constraint_block_rhs!(backend, block, 0),
        () -> update_constraint_variable_coefficients!(backend, 1, constraint),
        () -> backend_update_objective_coefficients!(backend, objective),
    )
    for operation in operations
        @test_throws ArgumentError operation()
    end
    mktempdir() do dir
        path = joinpath(dir, "unsupported.lp")
        @test_throws ArgumentError backend_write_to_file(backend, path)
        @test !ispath(path)
    end
    @test backend.status === :not_solved
    @test isnan(backend.objective_value)
    @test variable.lower_bounds == [0.0]
    @test variable.upper_bounds == [1.0]
    @test constraint.rhs == 3.0
end

@testset "CPLEX termination semantics (license-free contract)" begin
    for (raw, expected) in (
        (1, :optimal), (2, :unbounded), (3, :infeasible),
        (4, :infeasible_or_unbounded), (10, :iteration_limit),
        (11, :time_limit), (25, :time_limit), (12, :objective_limit),
        (13, :interrupted), (6, :numerical_error), (23, :feasible),
        (0, :other_error), (999, :other_error),
    )
        @test RM._cplex_status(false, raw) === expected
    end
    for (raw, expected) in (
        (101, :optimal), (102, :optimal), (103, :infeasible),
        (118, :unbounded), (119, :infeasible_or_unbounded),
        (133, :infeasible_or_unbounded), (107, :time_limit),
        (108, :time_limit), (131, :time_limit), (132, :time_limit),
        (105, :node_limit), (106, :node_limit), (104, :solution_limit),
        (111, :memory_limit), (113, :interrupted), (127, :feasible),
        (0, :other_error), (999, :other_error),
    )
        @test RM._cplex_status(true, raw) === expected
    end
end

@testset "Xpress termination semantics (license-free contract)" begin
    for (is_mip, raw, stop, expected) in (
        (false, 1, 0, :optimal), (false, 2, 0, :infeasible),
        (false, 5, 0, :unbounded), (false, 3, 0, :objective_limit),
        (false, 7, 0, :numerical_error), (false, 4, 0, :other_error),
        (false, 0, 0, :other_error), (false, 999, 0, :other_error),
        (true, 6, 0, :optimal), (true, 5, 0, :infeasible),
        (true, 4, 0, :feasible), (true, 7, 0, :infeasible_or_unbounded),
        (true, 3, 0, :other_error), (true, 999, 0, :other_error),
        (true, 6, 5, :optimal), (true, 6, 10, :optimal),
        (true, 4, 1, :time_limit), (true, 3, 1, :time_limit),
        (false, 4, 1, :time_limit), (false, 4, 2, :interrupted),
        (true, 4, 3, :node_limit), (false, 4, 4, :iteration_limit),
        (true, 4, 6, :solution_limit), (true, 4, 8, :memory_limit),
        (true, 4, 9, :interrupted), (false, 4, 13, :numerical_error),
        (true, 4, 14, :work_limit), (false, 1, 999, :other_error),
    )
        @test RM._xpress_status(is_mip, raw, stop) === expected
    end
end

end # module CommercialBackendContractTests
