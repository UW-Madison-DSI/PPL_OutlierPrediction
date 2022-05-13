using  Gen
using Plots

include("../step01-importing-data/utilities/read-files.jl")
DF = ReadDF("../../data/processed/DetrendedCov.csv")
xs = DF.Date
ys = DF.DetrendedN1
include("utilities/visualize.jl")


```
Creates a random linear model with some probablity any given point is an outlier. The models have the from:
    ys = slope*xs + intercept with added noise and some points being an outlier with no mean and large amount of noise

#Arguments
- 'xs::Vector{<:Real}' input data to be used to generate output.
```
@gen function Linear_regression_with_outliers(xs::Vector{<:Real})
    # First, generate the core parameters of the model. 
    #These parameters are random because we will MCMC to find what parameters have higher probability
    slope ~ normal(0, 1000)# We use normal because we want a symetric distribution around zero that has all values have positive probability
    intercept ~ normal(0, 1000)#Same as why we use normal for slope
    noise ~ gamma(250, 250)#noise needs to be always positive and we want it to look vaguely normal so gamma makes sense
    prob_outlier ~ uniform(0, 1)#no bias in the % of outlier
    
    #Used to set a 1 to 1 correspondence to ys
    n = length(xs)
    # Next, we generate the container for the y values.
    ys = Float64[]
    
    for i = 1:n
        # Decide whether this point is an outlier, and set
        # mean and standard deviation accordingly
        if ({:data => i => :is_outlier} ~ bernoulli(prob_outlier))
            (mu, std) = (0., 1e6)
        else
            (mu, std) = (xs[i] * slope + intercept, noise)
        end
        # Sample a y value for this point
        push!(ys, {:data => i => :y} ~ normal(mu, std))
    end
    ys
end;


function serialize_trace(trace)
    (xs,) = Gen.get_args(trace)
    n = length(xs)
    Dict(
         :points => zip(xs, [trace[:data => i => :y] for i in 1:n]),
         :outliers => [trace[:data => i => :is_outlier] for i in 1:n],
         :xs => xs,
         :ys => [xs[i] * trace[:slope] + trace[:intercept] for i in 1:n])
end



function visualize_trace(trace::Dict; title="")

    #Graph points
    outliers = [pt for (pt, outlier) in zip(trace[:points], trace[:outliers]) if outlier]
    inliers =  [pt for (pt, outlier) in zip(trace[:points], trace[:outliers]) if !outlier]
    scatter(map(first, inliers), map(last, inliers), markercolor="blue", label=nothing) 
    scatter!(map(first, outliers), map(last, outliers), markercolor="red", label=nothing)

    #Graph Line
    PLT = plot!(trace[:xs], trace[:ys], color = "black", lw = 3, label = nothing)


    #LogPlot
    scatter(map(first, inliers), log_modulus.(map(last, inliers)), markercolor="blue", label=nothing) 
    scatter!(map(first, outliers), log_modulus.(map(last, outliers)), markercolor="red", label=nothing)
    ExpPLT = plot!(trace[:xs], log_modulus.(trace[:ys]), color = "black", lw = 3, label = nothing)
    
    DuoPlot = plot(PLT,ExpPLT,title=title)
    return DuoPlot
end


function VizGenModel(GenFunction::DynamicDSLFunction)
    ts = collect(range(1, stop=400, length=400))#create equaly spaced input values, simmilar in shape to the ones we use
    traces = simulate(GenFunction, (ts,));#Creating model with no constraints on parameters
    visualize_trace(serialize_trace(traces))
end

show(VizGenModel(Linear_regression_with_outliers),"step02_test")


function make_constraints(ys::Vector{Float64})
    constraints = Gen.choicemap()
    for i=1:length(ys)
        constraints[:data => i => :y] = ys[i]
    end
    constraints
end;

observations = make_constraints(ys);


# Perform a single block resimulation update of a trace.
function block_resimulation_update(tr)
    # Block 1: Update the line's parameters
    line_params = select(:noise, :slope, :intercept)
    (tr, _) = mh(tr, line_params)
    
    # Blocks 2-N+1: Update the outlier classifications
    (xs,) = get_args(tr)
    n = length(xs)
    for i=1:n
        (tr, _) = mh(tr, select(:data => i => :is_outlier))
    end
    
    # Block N+2: Update the prob_outlier parameter
    (tr, _) = mh(tr, select(:prob_outlier))
    
    # Return the updated trace
    tr
end;

function VizGenMCMC(GenFunction::DynamicDSLFunction,xs,observations::DynamicChoiceMap,updateTrace::Function, NumFrame::Int64)
    t, = generate(GenFunction, (xs,), observations)#Create initial set of parameters to iterate on
    viz = @animate for i in 1:NumFrame
        t = updateTrace(t)
        visualize_trace(serialize_trace(t); title="Iteration $i/$NumFrame")
    end;
    gif(viz)
end
show(VizGenMCMC(Linear_regression_with_outliers, xs, observations,block_resimulation_update,100,RetAni=true),"step03.gif")