"""
    Returns a symetrical log that works for positive or negitive values

    # Arguments
    - 'x' - real number
    - 'm' - fudge factor that controls how much space the small values have
"""
function log_modulus(x,m = 10000)
    return sign(x)*(log(abs(x)+m)-log(m))
end


using Plots: scatter, scatter!, plot, plot!
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
    scatter(map(first, inliers), map(last, inliers), markercolor="blue", title ="Regular", label=nothing) 
    scatter!(map(first, outliers), map(last, outliers), markercolor="red", label=nothing)

    #Graph Line
    PLT = plot!(trace[:xs], trace[:ys], color = "black", lw = 3, label = nothing)


    #LogPlot
    scatter(map(first, inliers), log_modulus.(map(last, inliers)), markercolor="blue", title ="Log", label=nothing) 
    scatter!(map(first, outliers), log_modulus.(map(last, outliers)), markercolor="red", label=nothing)
    ExpPLT = plot!(trace[:xs], log_modulus.(trace[:ys]), color = "black", lw = 3, label = nothing)
    
    DuoPlot = plot(PLT,ExpPLT, plot_title = title)
    return DuoPlot
end



"""
    Creates a plot of the model creating a random scenario.

# Arguments
- 'GenFunction::DynamicDSLFunction' - @Gen model function that creates a random kind of output
"""
function VizGenModel(GenFunction::DynamicDSLFunction)
    ts = collect(range(1, stop=400, length=400))#create equaly spaced input values, simmilar in shape to the ones we use
    traces = simulate(GenFunction, (ts,));#Creating model with no constraints on parameters
    visualize_trace(serialize_trace(traces))
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
    Creates a GIF of the model MCMC permutation on the data

# Arguments
- 'GenFunction::DynamicDSLFunction' - @Gen model function with random parameters
- 'xs,observations::DynamicChoiceMap' - contstraints on the model. Forces the output to equal the real messurements
- 'updateTrace::Function' - A function oh how to MCMC chunk update the code
- 'NumFrame::Int64' - The number of frames to render
- 'RetEnd' - used to return the animation object instead of the gif
"""
function VizGenMCMC(GenFunction::DynamicDSLFunction,xs,observations::DynamicChoiceMap,updateTrace::Function, NumFrame::Int64; RetEnd::Bool=false)
    t, = generate(GenFunction, (xs,), observations)#Create initial set of parameters to iterate on
    viz = @animate for i in 1:NumFrame
        t = updateTrace(t)
        visualize_trace(serialize_trace(t); title="Iteration $i/$NumFrame")
    end;
    return(gif(viz))
end
