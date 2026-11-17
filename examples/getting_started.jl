using ReactiveModel

@model Dispatch(data, point) begin
    @variable production[g in 1:2] in [0.0, data.capacity[g]]
    @constraint balance sum(production[g] for g in 1:2) == data.demand[point.t]
    @objective Min sum(data.cost[g] * production[g] for g in 1:2)
end

data = (capacity=[80.0, 60.0], cost=[40.0, 70.0], demand=[100.0, 110.0])
problem = LinearProblem(ConstData(Dispatch(data), data), (t=1,);
                        backend=HiGHSBackend(), safe_mode=true)

first_cost = objective_value(problem, (t=1,))
solution = solve(problem, (t=2,))
second_cost = objective_value(solution)
production = [primal(solution, :production, g) for g in 1:2]

@assert isapprox(first_cost, 4600.0)
@assert isapprox(second_cost, 5300.0)
@assert production ≈ [80.0, 30.0]
(first_cost, second_cost, production)
