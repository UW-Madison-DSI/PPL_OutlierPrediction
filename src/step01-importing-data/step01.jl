using Plots, CSV, Dates
using DataFrames: DataFrame

#read data files using DataArrays
dataframe = CSV.read("../../data/Proccesed/DetrendedCov.csv", DataFrame)
dataframe[!,:Date] = convert.(Base.Float64,dataframe[!,:Date])

# display plot
xs = dataframe."Date"
ys = dataframe."DetrendedN1"
display(Plots.scatter(xs, ys, color="black", xlabel="X", ylabel="Y", label=nothing, title="Observations"))

# wait for keyboard input
readline()