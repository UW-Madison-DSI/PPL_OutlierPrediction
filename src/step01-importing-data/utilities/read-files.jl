using CSV: read
using DataFrames: DataFrame

"""
This function reads covid wastewater data from a file into a dataframe.
It coverts the N1 column to floats and converts the Date column to dates.
Other columns are read as strings.

# Arguments
- 'filePath::String' - the path to the file to read data from.
"""
function ReadDF(filePath::String)

    # use CSV package to load the CSV file
    dataframe =read(filePath, DataFrame)

    # convert N1 column to floats
    dataframe[!,:N1] = convert.(Base.Float64,dataframe[!,:N1])

    # convert Date column to dates
    dataframe[!,:Date] = convert.(Base.Float64,dataframe[!,:Date])

    return(dataframe)
end