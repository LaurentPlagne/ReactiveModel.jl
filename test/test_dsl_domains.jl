using Test
using ReactiveModel

@model DSLDomainBinary(data, point) begin
    @variable active[i in 1:data.n] in [0.0, 1.0] domain=:binary
    @constraint minimum sum(active[i] for i in 1:data.n) >= point.minimum
    @objective Min sum(active[i] for i in 1:data.n)
end

@model DSLDomainInteger(data, point) begin
    @variable units in [0.0, data.limit] domain=:integer
    @constraint minimum units >= point.minimum
    @objective Min units
end

@testset "DSL integer and binary domains" begin
    binary_data = (n=3,)
    definition = DSLDomainBinary(binary_data)
    @test definition.variables.active.domain === :binary
    lp = LinearProblem(ConstData(definition, binary_data), (minimum=1.5,))
    @test objective_value(lp, (minimum=1.5,)) ≈ 2.0
    sol = solve(lp, (minimum=2.5,))
    @test objective_value(sol) ≈ 3.0
    @test all(i -> primal(sol, :active, i) == 1.0, 1:3)

    integer_data = (limit=10.0,)
    integer_definition = DSLDomainInteger(integer_data)
    @test integer_definition.variables.units.domain === :integer
    integer_lp = LinearProblem(ConstData(integer_definition, integer_data), (minimum=2.2,))
    @test objective_value(integer_lp, (minimum=2.2,)) ≈ 3.0
    @test objective_value(integer_lp, (minimum=4.2,)) ≈ 5.0

    for (def, data) in ((definition, binary_data), (integer_definition, integer_data))
        clp = LinearProblem(ConstData(def, data), (minimum=1.5,); backend=ClpBackend())
        @test_throws ArgumentError solve(clp, (minimum=1.5,))
    end
    @test_throws ArgumentError ReactiveModel._parse_variable_domain([:(domain=:invalid)])
    @test_throws ArgumentError ReactiveModel._parse_variable_domain([:(domain=:binary), :(domain=:integer)])
end
