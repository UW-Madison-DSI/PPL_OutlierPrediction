using  Plots, CSV
using DataFrames: DataFrame

#using CSV package to load the CSV file
dataframe = CSV.read("../../Data/Proccesed/DetrendedCov.csv", DataFrame)
first(dataframe, 5)

#Plot the signal
Plots.scatter(dataframe."Date", dataframe."N1", color="black", xlabel="X", ylabel="Y", 
              label=nothing, title="Observations - regular data and outliers")

#Plot the detrended signal
Plots.scatter(dataframe."Date", dataframe."DetrendedN1", color="black", xlabel="X", ylabel="Y", 
              label=nothing, title="Observations - detrended data with outliers,")