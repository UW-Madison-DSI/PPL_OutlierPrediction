using  Gen
using Plots

include("../step01-importing-data/utilities/display.jl")
include("../step01-importing-data/utilities/read-files.jl")
DF = ReadDF("../../data/processed/DetrendedCov.csv")
xs = DF.Date
ys = DF.DetrendedN1
include("utilities/visualize.jl")


"""
Creates a random linear model with some probablity any given point is an outlier. The models is linear with added noise
and some points being an outlier with no mean and large amount of noise 

# Arguments
- 'xs::Vector{<:Real}' input data to be used to generate output.
"""
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




"""
Extract the infomation needed to plot from the more complex Gen trace object

#Arguments
- 'trace::Gen.DynamicDSLTrace' Gen trace of infomation about the model
"""
function serialize_trace(trace::Gen.DynamicDSLTrace)
    (xs,) = Gen.get_args(trace)
    n = length(xs)
    Dict(
         :points => zip(xs, [trace[:data => i => :y] for i in 1:n]),
         :outliers => [trace[:data => i => :is_outlier] for i in 1:n],
         :xs => xs,
         :ys => [xs[i] * trace[:slope] + trace[:intercept] for i in 1:n])
end



"""
    Reads a gen trace and outputs a plot showing the data with and without a log modulus transformation

# Arguments
- 'trace::Dict' - Dict containing the x values, the point y values, whether the data was flaged an outlier, and the model y values
- 'title::String' - Title of plot
"""
function visualize_trace(trace::Dict; title::String="")

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
    
    DuoPlot = plot(PLT,ExpPLT, plot_title = title)
    return DuoPlot
end


"""
    Creates a choicemap that forces the output of the model to be the true output of the data. 
    This is used to find the conditional probability that any model leads to this outcome

# Arguments
- 'ys::Vector{Float64}' - The output of the signal we are messuring
"""
function make_constraints(ys::Vector{Float64})
    constraints = Gen.choicemap()
    for i=1:length(ys)
        constraints[:data => i => :y] = ys[i]
    end
    constraints
end;








"""
    Perform a MCMC update of the Gen model updating. updates the global parameters the the local ones

# Arguments
- 'tr::Gen.DynamicDSLTrace' - The model trace containing the parameters that we update
"""
function block_resimulation_update(tr::Gen.DynamicDSLTrace)
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

"""
    Creates a GIF of the model MCMC permutation on the data

# Arguments
- 'GenFunction::DynamicDSLFunction' - @Gen model function with random parameters
- 'xs,observations::DynamicChoiceMap' - contstraints on the model. Forces the output to equal the real messurements
- 'updateTrace::Function' - A function oh how to MCMC chunk update the code
- 'NumFrame::Int64' - The number of frames to render
- 'RetAni' - used to return the animation object instead of the gif
"""
function VizGenMCMC(GenFunction::DynamicDSLFunction,xs,observations::DynamicChoiceMap,updateTrace::Function, NumFrame::Int64; RetAni::Bool=false)
    t, = generate(GenFunction, (xs,), observations)#Create initial set of parameters to iterate on
    viz = @animate for i in 1:NumFrame
        t = updateTrace(t)
        visualize_trace(serialize_trace(t); title="Iteration $i/$NumFrame")
    end;
    if(RetAni)
        return(viz)
    else
        return(gif(viz))
    end
end


#Main
show(VizGenModel(Linear_regression_with_outliers),"step02_test")
observations = make_constraints(ys);
show(VizGenMCMC(Linear_regression_with_outliers, xs, observations,block_resimulation_update,100,RetAni=true),"step03.gif")