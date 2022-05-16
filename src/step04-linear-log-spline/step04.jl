#===============================================================================
|                                                                              |
|                                   step04.jl                                  |
|                                                                              |
|==============================================================================|
|                                                                              |
|        This step generates and displays a linear segmented model with        |
|        the data mapped to a log scale.                                       |
|                                                                              |
|        Author(s): Steve Goldstein, Marlin Lee, Wansoo Cho,                   |
|                   AbeMegahed                                                 |
|                                                                              |
|        This file is subject to the terms and conditions defined in           |
|        'LICENSE.txt', which is part of this source code distribution.        |
|                                                                              |
|==============================================================================|
|     Copyright (C) 2022, Data Science Institute, University of Wisconsin      |
|==============================================================================#

using Gen
using Plots

include("../step01-importing-data/utilities/display.jl")
include("../step01-importing-data/utilities/read-files.jl")
include("../step02-linear-model/utilities/visualize.jl")

DF = ReadDF("../../data/processed/DetrendedCov.csv")
xs = DF.Date
ys = DF.N1
SubChunkSize = 20


"""
    Calculate what chunk point i is in.

    # Arguments
- 'i::Int' - index of input we want to find the associated chunk of
-  `chunkSize` - The size of the chunks.
"""
function DiffrenceIndex(i::Int; chunkSize = SubChunkSize)#helper function to find out what chunk a point is in
    return(div(i,chunkSize,RoundUp))
end


"""
    Calculate the y value the model generates before noise or outliers are added. Each chunk has its own slope and has a y intercept defined such that each chunk edge matchs. exact same function as in step 3.

    # Arguments
- `xs::Vector{Float64}` - The input vector that we are generating based on
- `Buffer_y::Float64` - The initial y intercept of that model
- `Slopes::Vector{Float64}` - The slopes of each chunk
"""
function yValCalc(xs::Vector{Float64}, Buffer_y::Float64, Slopes::Vector{Float64})
    n = length(xs)
    NumChunks = DiffrenceIndex(n)


    # Calculating the 'y intercept' of each chunk to make sure each line connects to the last one.
    # Because each intercept gets added to the last one we take the cumalitive sum to get the total ofset needed at each step.
    # The first value should be the initial ofset Buffer_y to get everything aligned.
    # ysOfseted = [Buffer_y, Slope[chunk](x[chunk]- x[Last chunk])]
    ysOfseted = cumsum(pushfirst!([Slopes[i]*(xs[(i)*SubChunkSize] - xs[(i-1)*SubChunkSize+1]) for i=1:(NumChunks-1)],Buffer_y))
    
    # Calculates the change of y from the previous chunk to the current x. We combine this with a set of y offset values in the next step to get the true mu fed into the normal distribution.
    # TrueDeltaMu n = Slope[chunk](x[i]- x[Last chunk])
    TrueDeltaMu = [Slopes[DiffrenceIndex(i)]*(xs[i] - xs[div(i-1,SubChunkSize,RoundDown)*SubChunkSize+1]) for i=1:n]
    ys = [TrueDeltaMu[i] + ysOfseted[DiffrenceIndex(i)] for i=1:n]
end


"""
    Creates a random model where the data is broken into chunks and a Linear line is fit on the data with noise and some probability that they are outliers.

    Same type of model generater as step 3 but with the distributions changed such that its modeling the log of the data instead.

    # Arguments
- `xs::Vector{Float64}` - The input vector that we are generating based on
"""
@gen function Log_Linear_Spline_with_outliers(xs::Vector{<:Real})
    # First we calculate some useful values needed for the list comprehension in the next steps.
    n = length(xs)
    NumChunks = DiffrenceIndex(n)

    # Next, we generate some parameters of the model. There are three types of randomly made perameters. First are the constant ones that are unique to the process. These are generated first.
    # Second are the ones that are unique to the individual chunks. These loop from 1 to NumChunks.
    # Last are the ones that vary for every point. These range from 1 to n.

    # Unique to process

    # Where the series starts. In the log model this is around 12 and I give it a pretty big window
    Buffer_y ~ normal(10, 6) 
    
    # The probability any given point is a outlier
    prob_outlier ~ uniform(.05, .1)

    # The scaling factor on outliers:
    OutlierDeg ~ uniform(2, 5)
    
    # Unique to chunk

    # The data apears to have no slope over 3 so a sd of 2 should capture the true slopes with high probability
    Slopes = [{(:slope, i)} ~ normal(0, .01) for i=1:NumChunks]

    # The distribution of the noise. It gets fed into the sd of a normal distribution so the distribution of the noise needs to be always positive
    noise = [{(:noise, i)} ~ gamma(2, .15) for i=1:NumChunks]
   
    # EveryPoint

    # Is using the prob_outlier vector above to decide if each point is an outlier. the model we are using now has the slope and sd $OutlierDeg times larger then the non outliers so we times the mu and sd by this value in the last step.
    PointOutlier = ((OutlierDeg-1)*[{:data => i => :is_outlier} ~ bernoulli(prob_outlier) for i=1:n] .+ 1)

    # The random var fit to the actual data. It is created as a combination of previous parts.
    # The process was discribed in previous steps.
    # ys = normal(mu, sd)
    TrueVec = yValCalc(xs,Buffer_y,Slopes)
    ys = [{:data => i => :y} ~ normal(
        TrueVec[i],
        noise[DiffrenceIndex(i)]*PointOutlier[i]
        ) for i=1:n]
    ys
end;


"""
    Extract the infomation needed to plot from the more complex Gen trace object. exp the y values so we view the data in the same form as the other sections.

    # Arguments
- 'trace::Gen.DynamicDSLTrace' Gen trace of infomation about the model.
"""
function serialize_trace(trace::Gen.DynamicDSLTrace)
    (xs,) = Gen.get_args(trace)
    n = length(xs)
    NumChunks = div(n, SubChunkSize, RoundUp)
    slopes = [trace[(:slope, i)] for i in 1:NumChunks]
    FlatDict = Dict(
          :points => zip(xs, [exp(trace[:data => i => :y]) for i in 1:n]),
          :outliers => [trace[:data => i => :is_outlier] for i in 1:n],
         :xs => xs,
         :ys => exp.(yValCalc(xs, trace[:Buffer_y], slopes)))
    return(FlatDict)
end


"""
    Perform a MCMC update of the Gen model updating.  Updates the global parameters the the local ones.

    # Arguments
- 'tr::Gen.DynamicDSLTrace' - The model trace containing the parameters that we update.
"""
function block_resimulation_update(tr::Gen.DynamicDSLTrace)
    (xs,) = get_args(tr)
    n = length(xs)
    NumChunks = div(n, SubChunkSize, RoundUp)

    for j=1:(NumChunks)

        # Block 1: Update the line's parameters.
        line_params = select((:noise,j-1), (:slope,j-1),
                    (:noise,j), (:slope,j),
                    (:noise,j+1), (:slope,j+1) )
        (tr, _) = mh(tr, line_params)
    end
    
    # Blocks 2-N+1: Update the outlier classifications.
    for i=1:n
        (tr, _) = mh(tr, select(:data => i => :is_outlier))
    end
    
    (tr, _) = mh(tr, select(:prob_outlier, :OutlierDeg, :Buffer_y))
    
    # Return the updated trace.
    tr
end;


# Main
show(VizGenModel(Log_Linear_Spline_with_outliers),"step04_test.png")

# Shows a gif of the MCMC working on the Waste Water data
observations = make_constraints(log.(ys));
show(VizGenMCMC(Log_Linear_Spline_with_outliers, xs, observations, block_resimulation_update, 100),"step05.gif")
