using Sockets;

"""
This function is used to determine if the current Julia program 
is running on the local computer or on a remote host.
"""
function is_local() 
    return startswith(string(Sockets.getipaddr()), "192.168")
end

"""
This function is used to show a plot either in a popop window
when running locally or saved to a file when running remotely.
"""
function show(plot, filename)
    if (is_local())

        # Render the Animation object into a GIF object
        if (typeof(plot) == Animation)
            gif(plot, filename)
            println("Saved plot to $filename.")
        else
            # Open a new window to display the plot.
            display(plot)

            # Wait for keyboard input to close plot window.
            readline()
        end
    else
        if (typeof(plot) == Animation)
            gif(plot, filename)
        else
            # Save plot to an image file.
            savefig(filename)
        end

        # Notify the current user that plot has been saved.
        println("Saved plot to $filename.")
    end
end

#=
Run in "headless" mode if remote.  This prevents various
plotting related error messages from being displayed when
running on a remote server.
=#
if (!is_local())
        println("Running in headless mode...");
        ENV["GKSwstype"] = "100"
end
