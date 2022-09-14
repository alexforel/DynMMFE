function demand_distribution_pattern(demandPattern::String, simulationLength::Integer)
    ## demand_distribution_pattern(demandPattern::String, simulationLength::Integer)
    # Create demand distributions following the user-specified pattern over the simulation length. Note that the Stationary and Seasonal patterns are stationary, even if specified as distributions because they have 0 variance.
    # Input:
    #   - demandPattern: String, user-specified demand pattern
    #   - simulationLength: Int, number of time periods in simulation
    # Output:
    #   - demandDistributions: Matrix[Distributions], list of distributions over the simulationLength periods

    demandDistributions = Matrix(undef, simulationLength, 1)

    # Define demand pattern
    if demandPattern == "Stationary"
        demandDistributions .= Normal(100, 0)
    elseif demandPattern == "Random"
        demandDistributions .= Uniform(50, 150)
    elseif demandPattern == "Seasonal"
        seasonLength = T # Demand is seasonal with periodicity simulationLength
        meanDem = zeros(seasonLength)
        # Seasonal demand generated from sinusoid function
        for t in 1:seasonLength
            meanDem[t] = 160 * sin(0.1 + pi * (t-1) / seasonLength)
        end
        for t in 1:simulationLength
            if t <= seasonLength
                demandDistributions[t] = Normal(meanDem[t], 0)
            else
                demandDistributions[t] = Normal(meanDem[t - seasonLength], 0)
            end
        end
    else
        error("Wrong demand pattern. User specified \"$demandPattern\" but only \"Stationary\", \"Random\", or \"Seasonal\" are supported.")
    end

    # Output
    return demandDistributions
end
