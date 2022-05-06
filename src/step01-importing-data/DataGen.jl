using  Plots, CSV
using DataFrames: DataFrame

function ReadDF(FilePath::String = "../../Data/Proccesed/DetrendedCov.csv")
    #using CSV package to load the CSV file
    dataframe = CSV.read("../../Data/Proccesed/DetrendedCov.csv", DataFrame)
    dataframe[!,:N1] = convert.(Base.Float64,dataframe[!,:N1])
    dataframe[!,:Date] = convert.(Base.Float64,dataframe[!,:Date])
    return(dataframe)
end

dataframe = ReadDF()
first(dataframe, 5)

#Plot the signal
Plots.scatter(dataframe."Date", dataframe."N1", color="black", xlabel="X", ylabel="Y", 
              label=nothing, title="Observations - regular data and outliers")

#Plot the detrended signal
Plots.scatter(dataframe."Date", dataframe."DetrendedN1", color="black", xlabel="X", ylabel="Y", 
              label=nothing, title="Observations - detrended data with outliers,")