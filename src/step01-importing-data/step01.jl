include("./utilities/read-files.jl")

# Read the data from a CSV file.
dataframe = ReadDF("../../Data/Proccesed/DetrendedCov.csv")
first(dataframe, 5)

# Plot the concentration data (N1 column) in red.
display(Plots.scatter(dataframe."Date", dataframe."N1", color="red", xlabel="Date", ylabel="Concentration", 
label="N1", title="Observations"))

# Plot the detrended concentration data (DetrendedN1 column) in blue.
display(Plots.scatter!(dataframe."Date", dataframe."DetrendedN1", color="blue", xlabel="Date", ylabel="Concentration", 
label="N1 (Detrended)", title="Observations"))

# Wait for keyboard input
readline()