#===============================================================================
|                                                                              |
|                                   step01.jl                                  |
|                                                                              |
|==============================================================================|
|                                                                              |
|        This step reads in and displays the covid wastewater data.            |
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

using Plots

include("./utilities/read-files.jl")
include("./utilities/display.jl")

# Read the data from a CSV file.
dataframe = ReadDF("../../data/processed/DetrendedCov.csv")
first(dataframe, 5)

# Plot the concentration data (N1 column) in red.
Plots.scatter(dataframe."Date", dataframe."N1", color="red", xlabel="Date", ylabel="Concentration", 
label="N1", title="Observations")

# Plot the detrended concentration data (DetrendedN1 column) in blue.
show(Plots.scatter!(dataframe."Date", dataframe."DetrendedN1", color="blue", xlabel="Date", ylabel="Concentration", label="N1 (Detrended)", title="Observations"), "step01.png")