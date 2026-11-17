# Standalone smoke test: run in a fresh Julia process/environment.
# Do not include this from a session which has already imported comparison tools.
using Test
using ReactiveModel

@testset "Standalone ReactiveModel uses HiGHS_jll without MOI" begin
    @test !any(id -> id.name in ("HiGHS", "MathOptInterface"), keys(Base.loaded_modules))
    @test any(lib -> lib.interface === :lp64,
        ReactiveModel._HiGHSCAPI.LinearAlgebra.BLAS.lbt_get_config().loaded_libs)
    definition = ProblemDefinition(
        variables=(x=VariableBlock(:x; lower=1.0, upper=2.0),),
        constraints=NamedTuple(),
        objective=Objective(:min, (eq, d, p) -> set_coeff!(eq, :x, 1.0)),
    )
    lp = LinearProblem(definition, nothing; backend=HiGHSBackend())
    @test objective_value(lp, nothing) ≈ 1.0
    @test !any(id -> id.name in ("HiGHS", "MathOptInterface"), keys(Base.loaded_modules))
end
