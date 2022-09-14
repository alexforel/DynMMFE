function average_and_relative_costs_per_simulation_environment(inputResultDf)
    ## average_and_relative_costs_per_simulation_environment(inputResultDf)
    # Calculate the average absolute and relative costs of the MMFE models w.r.t. the deterministic model. Also evaluates the statistical significance of the relative cost being different than 100\%.
    # Input:
    #   - inputResultDf
    # Output:
    #   - averageCost, relativeCost, statSignRelativeCost

    # Initialise matrices
    averageCost = zeros(length(demandPatternVector) * length(uncertaintyCoeffVector) * length(capacityVector), NB_MODELS)
    relativeCost = zeros(length(demandPatternVector) * length(uncertaintyCoeffVector) * length(capacityVector), NB_MODELS)
    statSignRelativeCost = zeros(length(demandPatternVector) * length(uncertaintyCoeffVector) * length(capacityVector) + 1, NB_MODELS-1)

    # For each parameter combination, filter the dataframe to them
    for demP in demandPatternVector
        tempDf = filter(row -> row[:demandPattern] == demP, inputResultDf)
        for u in uncertaintyCoeffVector
            tempDf2 = filter(row -> row[:uncertaintyCoeff] == u, tempDf)
            for c in capacityVector
                tempDf3 = filter(row -> row[:capacity] == c, tempDf2)
                # Get absolute costs of each solution method
                detCost = filter(row -> row[:model] == "deterministic", tempDf3)[!, :cost]
                ddstochCost = filter(row -> row[:model] == "dd-stoch", tempDf3)[!, :cost]
                addCost = filter(row -> row[:model] == "additive", tempDf3)[!, :cost]
                multCost = filter(row -> row[:model] == "multiplicative", tempDf3)[!, :cost]
                extendedCost = filter(row -> row[:model] == "extended", tempDf3)[!, :cost]
                # Store mean
                a = findfirst(isequal(demP), demandPatternVector)
                b = findfirst(isequal(u), uncertaintyCoeffVector)
                c = findfirst(isequal(c), capacityVector)
                rowToStore = length(capacityVector) * length(uncertaintyCoeffVector) * (a-1) + length(capacityVector) * (b-1) + c
                meanCostVector = mean([detCost ddstochCost addCost multCost extendedCost], dims = 1)
                averageCost[rowToStore, :] = meanCostVector
                # Relative cost
                reldetCost = detCost ./ detCost * 100
                relddstochCost = ddstochCost ./ detCost * 100
                reladdCost = addCost ./ detCost* 100
                relmultCost = multCost ./ detCost* 100
                relextendedCost = extendedCost ./ detCost* 100
                relativeCost[rowToStore, :] = mean([reldetCost relddstochCost reladdCost relmultCost relextendedCost], dims = 1)
                # Assess statistical difference of improvement relative to deterministic
                statSignRelativeCost[rowToStore, 1] =  pvalue(OneSampleTTest(relddstochCost, 100)) <= 0.05
                statSignRelativeCost[rowToStore, 2] =  pvalue(OneSampleTTest(reladdCost, 100)) <= 0.05
                statSignRelativeCost[rowToStore, 3] =  pvalue(OneSampleTTest(relmultCost, 100)) <= 0.05
                statSignRelativeCost[rowToStore, 4] =  pvalue(OneSampleTTest(relextendedCost, 100)) <= 0.05
            end
        end
    end
    # Statistical significance over all instances
    detCost = filter(row -> row[:model] == "deterministic", inputResultDf)[!, :cost]
    ddstochCost = filter(row -> row[:model] == "dd-stoch", inputResultDf)[!, :cost]
    addCost = filter(row -> row[:model] == "additive", inputResultDf)[!, :cost]
    multCost = filter(row -> row[:model] == "multiplicative", inputResultDf)[!, :cost]
    extendedCost = filter(row -> row[:model] == "extended", inputResultDf)[!, :cost]
    reldetCost = detCost ./ detCost * 100
    relddstochCost = ddstochCost ./ detCost * 100
    reladdCost = addCost ./ detCost * 100
    relmultCost = multCost ./ detCost * 100
    relextendedCost = extendedCost ./ detCost * 100
    # Assess statistical difference of improvement relative to deterministic
    statSignRelativeCost[end, 1] =  pvalue(OneSampleTTest(relddstochCost, 100)) <= 0.05
    statSignRelativeCost[end, 2] =  pvalue(OneSampleTTest(reladdCost, 100)) <= 0.05
    statSignRelativeCost[end, 3] =  pvalue(OneSampleTTest(relmultCost, 100)) <= 0.05
    statSignRelativeCost[end, 4] =  pvalue(OneSampleTTest(relextendedCost, 100)) <= 0.05

    return averageCost, relativeCost, statSignRelativeCost
end
