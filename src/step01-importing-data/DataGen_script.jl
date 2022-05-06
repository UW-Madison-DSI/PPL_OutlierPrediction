function ReadDF(FilePath::String = "../../Data/Proccesed/DetrendedCov.csv")
    #using CSV package to load the CSV file
    dataframe = CSV.read("../../Data/Proccesed/DetrendedCov.csv", DataFrame)
    dataframe[!,:N1] = convert.(Base.Float64,dataframe[!,:N1])
    dataframe[!,:Date] = convert.(Base.Float64,dataframe[!,:Date])
    return(dataframe)
end