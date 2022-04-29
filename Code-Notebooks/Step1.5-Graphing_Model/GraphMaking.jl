using  Plots

function visualize_trace(trace::Trace; title="")
    trace = serialize_trace(trace)
    #Graph points
    outliers = [pt for (pt, outlier) in zip(trace[:points], trace[:outliers]) if outlier]
    inliers =  [pt for (pt, outlier) in zip(trace[:points], trace[:outliers]) if !outlier]
    PLT = Plots.scatter(map(first, inliers), map(last, inliers), markercolor="blue", label=nothing, title=title) 
    PLT = Plots.scatter!(map(first, outliers), map(last, outliers), markercolor="red", label=nothing)

    #Graph Line
    PLT = Plots.plot!(trace[:xs], trace[:ys], color = "black", lw = 3, label = nothing)
    return PLT
end