"""
This function is used to determine if the current Julia program 
is running on the local computer or on a remote host.
"""
function is_local() 
    return endswith(gethostname(), ".local")
end

"""
This function is used to show a plot either in a popop window
when running locally or saved to a file when running remotely.
"""
function show(plot, filename)
    if (is_local()) 
        display(plot)

        # Wait for keyboard input to close plot window.
        readline()
    else
        # Save plot to an image file.
        savefig(filename)

        # Notify the current user that plot has been saved.
        println("Saved plot to $filename.")
    end
end