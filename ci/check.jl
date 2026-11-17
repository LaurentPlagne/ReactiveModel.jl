# Run package or documentation checks from a standalone source copy. Resolution
# uses the selected Julia version without modifying the checkout's manifests.
#
#   julia --startup-file=no ci/check.jl test
#   julia --startup-file=no ci/check.jl docs
function check(mode::String)
    mode in ("test", "docs") || throw(ArgumentError("expected test or docs"))
    root = dirname(@__DIR__)
    stage = mktempdir(; prefix="reactivemodel-ci-")
    for name in ("Project.toml", "LICENSE", "src", "ext", "test", "examples")
        cp(joinpath(root, name), joinpath(stage, name))
    end
    julia = Base.julia_cmd()
    println("Julia ", VERSION, "; standalone source: ", stage)
    if mode == "test"
        code = "using Pkg; Pkg.instantiate(); Pkg.test(; allow_reresolve=false)"
        run(`$julia --startup-file=no --project=$stage -e $code`)
    else
        # Documentation excerpts read these files as text; benchmark packages
        # and their manifest are not needed to build the manual.
        mkpath(joinpath(stage, "benchmarks"))
        for name in ("valley_comparison.jl", "jump_user_model.jl")
            cp(joinpath(root, "benchmarks", name), joinpath(stage, "benchmarks", name))
        end
        docs = joinpath(stage, "docs")
        mkpath(docs)
        for name in ("Project.toml", "make.jl", "src")
            cp(joinpath(root, "docs", name), joinpath(docs, name))
        end
        run(`$julia --startup-file=no --project=$docs -e 'using Pkg; Pkg.instantiate()'`)
        run(`$julia --startup-file=no --project=$docs $(joinpath(docs, "make.jl"))`)
        # The checker locates docs/build relative to its own source file.
        mkpath(joinpath(stage, "audit"))
        checker = joinpath(stage, "audit", "check-documentation-links.py")
        cp(joinpath(root, "audit", "check-documentation-links.py"), checker)
        run(`python3 $checker`)
    end
    return nothing
end

length(ARGS) == 1 || throw(ArgumentError("usage: ci/check.jl test|docs"))
check(only(ARGS))
