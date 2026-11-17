using Test
using ReactiveModel

@testset "Optional backend availability and unsupported ClpJL operations" begin
    if Base.get_extension(ReactiveModel, :ReactiveModelTinyHiGHSExt) === nothing
        @test_throws ArgumentError TinyHiGHSBackend()
    end
    # These guards do not require the optional solver to be installed or initialized.
    backend = ClpJLBackend(nothing, nothing, nothing, false,
                          Float64[], Float64[], Float64[], Float64[], NaN, :not_solved)
    row = ReactiveModel._BackendConstraint([1], [2.0], 3.0, :le)
    @test_throws ArgumentError update_constraint_variable_coefficients!(backend, 1, row)
    mktempdir() do dir
        file = joinpath(dir, "unsupported.lp")
        @test_throws ArgumentError backend_write_to_file(backend, file)
        @test !ispath(file)
    end
end

# Run this same file from an environment with TinyHiGHS to exercise the extension.
if Base.find_package("TinyHiGHS") !== nothing
    import TinyHiGHS
    include("test_tinyhighs.jl")
else
    @testset "Optional TinyHiGHS integration (not installed)" begin
        @test_skip false
    end
end
