include("workload.jl")

function (@main)(args::Vector{String})
    n = isempty(args) ? 4 : parse(Int, only(args))
    return n > 0 && ReactiveModelTrimWorkload.verify(n) ? 0 : 1
end
