function visualize_trace(trace::Trace; title="", YscaleFunc::Function =  identity)
    trace = serialize_trace(trace)
    #Graph points
    outliers = [pt for (pt, outlier) in zip(trace[:points], trace[:outliers]) if outlier]
    inliers =  [pt for (pt, outlier) in zip(trace[:points], trace[:outliers]) if !outlier]
    PLT = Plots.scatter(map(first, inliers), YscaleFunc.(map(last, inliers)), markercolor="blue", label=nothing, title=title) 
    PLT = Plots.scatter!(map(first, outliers), YscaleFunc.(map(last, outliers)), markercolor="red", label=nothing)
    #Graph Line
    PLT = Plots.plot!(trace[:xs], YscaleFunc.(trace[:ys]), color = "black", lw = 3, label = nothing)
    return PLT
end

function make_constraints(ys::Vector{Float64})
    constraints = Gen.choicemap()
    for i=1:length(ys)
        constraints[:data => i => :y] = ys[i]
    end
    constraints
end;
