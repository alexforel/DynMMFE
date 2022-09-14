function print_message_iteration(simIndex)
    ## print_message_iteration(simIndex)
    # Print a dot every 10 simIndex and break line every 100.
    # Input:
    #   - simIndex
    # Output:

    if (simIndex % 10) == 0
        if (simIndex % 100) == 0
            Core.print(".\n")
        else
            Core.print(".")
        end
    end

    return nothing
end
