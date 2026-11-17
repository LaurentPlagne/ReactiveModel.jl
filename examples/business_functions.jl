module BusinessFunctionsExample
using ReactiveModel

# Deterministic business functions: both inputs are explicit. Their internal
# branches and calculations remain ordinary Julia, outside the modeling DSL.
function available_capacity(study, point)
    thermal_factor = clamp(1.0-study.derating*max(0.0,point.temperature-25.0),0.0,1.0)
    study.capacity * clamp(point.availability,0.0,1.0) * thermal_factor
end
function conversion_efficiency(study, point)
    temperature_factor = 1.0 / (1.0+0.001*(point.temperature-20.0)^2)
    study.efficiency * temperature_factor * (point.boost ? 1.1 : 1.0)
end
function requested_output(study, point)
    cooling = study.cooling*max(0.0,point.temperature-25.0)
    max(study.minimum_load,study.base_load*point.load_factor+cooling)
end
operating_cost(study,point) = study.fuel_price*point.price
fixed_cost(study,point) = point.boost ? study.boost_cost : 0.0

@model BusinessDispatch(study,point) begin
    @variable input in [0.0,available_capacity(study,point)]
    @constraint delivery conversion_efficiency(study,point)*input >= requested_output(study,point)
    @objective Min operating_cost(study,point)*input + fixed_cost(study,point)
end

# A deliberately incomplete condition, used to demonstrate a missed update.
@model IncompleteCapacity(study,point) begin
    @variable input in [0.0,available_capacity(study,point)] update_when :availability
    @objective Max input
end

study_data() = (capacity=30.0,derating=0.01,efficiency=0.9,cooling=0.1,
                minimum_load=1.0,base_load=10.0,fuel_price=2.0,boost_cost=3.0)
initial_point() = (temperature=20.0,availability=1.0,load_factor=0.7,boost=false,price=1.0)

function run_study(;backend=HiGHSBackend())
    data=study_data()
    problem=LinearProblem(ConstData(BusinessDispatch(data),data),initial_point();backend,safe_mode=true)
    values=Float64[]
    # A driver may branch, skip candidates or revisit points. No numerical
    # setters, model indices or update order appear in this control flow.
    for temperature in (20.0,35.0), load_factor in (0.0,0.7,1.1)
        load_factor==0.0 && continue
        for boost in (false,true)
            point=(temperature=temperature,availability=0.9,load_factor=load_factor,
                   boost=boost,price=1.2)
            value=objective_value(problem,point)
            required=requested_output(data,point)/conversion_efficiency(data,point)
            @assert required <= available_capacity(data,point)
            @assert isapprox(value,operating_cost(data,point)*required+fixed_cost(data,point))
            push!(values,value)
        end
    end
    @assert length(values)==8
    @assert all(c.missed==0 for c in update_report(problem))
    @assert any(c.unnecessary>0 for c in update_report(problem))
    (;problem,values)
end

function demonstrate_missed_update(;backend=HiGHSBackend())
    data=study_data()
    first=initial_point()
    hot=merge(first,(temperature=35.0,))
    problem=LinearProblem(ConstData(IncompleteCapacity(data),data),first;backend,safe_mode=true)
    @assert objective_value(problem,first)==30.0
    # Availability is unchanged. The predicate misses the temperature effect;
    # safe mode still applies the capacity change before solving.
    @assert objective_value(problem,hot)==27.0
    @assert only(c for c in update_report(problem) if c.component===:bounds).missed==1
    problem
end
end

if abspath(PROGRAM_FILE)==@__FILE__
    result=BusinessFunctionsExample.run_study()
    println("Checked objectives: ",result.values)
    foreach(println,BusinessFunctionsExample.update_report(result.problem))
    corrected=BusinessFunctionsExample.demonstrate_missed_update()
    println("Corrected incomplete condition:")
    foreach(println,BusinessFunctionsExample.update_report(corrected))
end
