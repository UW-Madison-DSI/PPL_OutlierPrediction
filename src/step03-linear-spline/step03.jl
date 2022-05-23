#===============================================================================
|                                                                              |
|                                   step03.jl                                  |
|                                                                              |
|==============================================================================|
|                                                                              |
|        This step generates and displays a segmented linear model.            |
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
SubChunkSize = 35   #Number of points per chunk


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
    Calculate the y value the model generates before noise or outliers are added. 

    Each chunk has its own slope and has a y intercept defined such that each chunk edge matchs.

    # Arguments
- `xs::Vector{Float64}` - The input vector that we are generating based on
- `Buffer_y::Float64` - The initial y intercept of that model
- `Slopes::Vector{Float64}` - The slopes of each chunk
"""
function yValCalc(xs::Vector{Float64}, Buffer_y::Float64, Slopes::Vector{Float64})
    n = length(xs)
    NumChunks = DiffrenceIndex(n)

    # Calculating the 'y intercept' of each chunk to make sure each line connects to the last one.
    # Because each intercept gets added to the last one we take the cumulative sum to get the total offset needed at each step.
    # The first value should be the initial offset Buffer_y to get everything aligned.
    # ysOfseted = [Buffer_y, Slope[chunk](x[chunk]- x[Last chunk])]
    ysOfseted = cumsum(pushfirst!([Slopes[i]*(xs[(i)*SubChunkSize] - xs[(i-1)*SubChunkSize+1]) for i=1:(NumChunks-1)],Buffer_y))
    
    # Calculates the change of y from the previous chunk to the current x. We combine this with a set of y offset values.
    # In the next step to get the true mu fed into the normal distribution
    # TrueDeltaMu n = Slope[chunk](x[i]- x[Last chunk])
    TrueDeltaMu = [Slopes[DiffrenceIndex(i)]*(xs[i] - xs[div(i-1,SubChunkSize,RoundDown)*SubChunkSize+1]) for i=1:n]
    ys = [TrueDeltaMu[i] + ysOfseted[DiffrenceIndex(i)] for i=1:n]
end


"""
    Creates a random model where the data is broken into chunks and a linear line is fit on the data with noise and some probability that they are outliers.

    # Arguments
- `xs::Vector{Float64}` - The input vector that we are generating based on
"""
@gen function Linear_Spline_with_outliers(xs::Vector{<:Real})
    # First we calculate some useful values needed for the list comprehension in the next steps.
    n = length(xs)
    NumChunks = DiffrenceIndex(n)

    # Next, we generate some parameters of the model. There are three types of randomly made perameters. 
    # First are the constant ones that are unique to the process. These are generated first.
    # Second are the ones that are unique to the individual chunks. These loop from 1 to NumChunks.
    # Last are the ones that vary for every point. These range from 1 to n.

    # Unique to process

    # Where the series starts. In the log model this is around 12 and I give it a pretty big window
    Buffer_y ~ gamma(200, 200)
    
    # The probability any given point is a outlier
    prob_outlier ~ uniform(.05, .1)

    # The scaling factor on outliers:
    OutlierDeg ~ uniform(1, 5)
    
    # Unique to chunk

    # The data apears to have no slope over 3 so a sd of 2 should capture the true slopes with high probability.
    Slopes = [{(:slope, i)} ~ normal(0, 3000) for i=1:NumChunks]

    # The distribution of the noise. It gets fed into the sd of a normal distribution so the distribution of the noise needs to be always positive.
    noise = [{(:noise, i)} ~ gamma(150, 150) for i=1:NumChunks]
   
    # EveryPoint is using the prob_outlier vector above to decide if each point is an outlier. the model we are using now has the slope and sd $OutlierDeg times larger then the non outliers so we times the mu and sd by this value in the last step.
    PointOutlier = ((OutlierDeg-1)*[{:data => i => :is_outlier} ~ bernoulli(prob_outlier) for i=1:n] .+ 1)

    TrueVec = yValCalc(xs,Buffer_y,Slopes)
    ys = [{:data => i => :y} ~ 
        normal(
            TrueVec[i] * PointOutlier[i],            #mean of normal rand var
            noise[DiffrenceIndex(i)]                 #var of normal rand var
        ) 
        for i=1:n]
    ys
end


"""
    Extract the infomation needed to plot from the more complex Gen trace object.

    # Arguments
- 'trace::Gen.DynamicDSLTrace' Gen trace of infomation about the model.
"""
# Get seralize trace to accept function instaed of unique code for each version.
function serialize_trace(trace::Gen.DynamicDSLTrace)
    (xs,) = Gen.get_args(trace)
    n = length(xs)
    NumChunks = div(n, SubChunkSize, RoundUp)
    slopes = [trace[(:slope, i)] for i in 1:NumChunks]
    FlatDict = Dict(
          :points => zip(xs, [trace[:data => i => :y] for i in 1:n]),
          :outliers => [trace[:data => i => :is_outlier] for i in 1:n],
         :xs => xs,
         :ys => yValCalc(xs, trace[:Buffer_y], slopes))
    return(FlatDict)
end


"""
    Perform a MCMC update of the Gen model updating. updates the global parameters the the local ones

    # Arguments
- 'tr::Gen.DynamicDSLTrace' - The model trace containing the parameters that we update
"""
function block_resimulation_update(tr::Gen.DynamicDSLTrace)
    (xs,) = get_args(tr)
    n = length(xs)
    NumChunks = div(n, SubChunkSize, RoundUp)
    for j=1:NumChunks
        # Block 1: Update the line's parameters.
        line_params = select((:noise,j), (:slope,j))
        (tr, _) = mh(tr, line_params)
    end
    
    # Blocks 2-N+1: Update the outlier classifications.
    for i=1:n
        (tr, _) = mh(tr, select(:data => i => :is_outlier))
    end
    
    # Block N+2: Update the prob_outlier parameter.
    (tr, _) = mh(tr, select(:prob_outlier, :Buffer_y, :OutlierDeg))
    
    # Return the updated trace.
    tr
end;

# Main

# Constrain the model so the output in the wastewater output.
observations = make_constraints(ys);

# Shows a gif of the MCMC working on the Waste Water data.
show(VizGenMCMC(Linear_Spline_with_outliers, xs, observations,block_resimulation_update,300),"step03.gif")
