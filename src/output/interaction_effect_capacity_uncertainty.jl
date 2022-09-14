function interaction_effect_capacity_uncertainty(inputString, inputDf)
    ## interaction_effect_capacity_uncertainty(inputDf)
    # Calculate the average relative cost of pla compared to det and its configuration interval for each simulation setting. Output the two measures in two matrices.
    # Input:
    #   - inputDf
    # Output:
    #   - averageRelativeCost, costConfInterval

    # Initialise matrices
    averageRelativeCost = zeros(length(uncertaintyCoeffVector), length(capacityVector)*length(demandPatternVector))
    costConfInterval = zeros(length(uncertaintyCoeffVector), length(capacityVector)*length(demandPatternVector))

    # For each simulation setting, filter to det and pla results. Calculate relative cost and confidence interval
    for d in 1:length(demandPatternVector)
        # Filter to demand pattern
        demFiltDf = filter(row -> row[:demandPattern] == demandPatternVector[d], inputDf)
        tempAverageRelativeCost = zeros(length(uncertaintyCoeffVector), length(capacityVector))
        tempCostConfInterval = zeros(length(uncertaintyCoeffVector), length(capacityVector))
        for u in 1:length(uncertaintyCoeffVector)
            tempDf = filter(row -> row[:uncertaintyCoeff] == uncertaintyCoeffVector[u], demFiltDf)
            for c in 1:length(capacityVector)
                tempDf2 = filter(row -> row[:capacity] == capacityVector[c], tempDf)
                filteredDetCost = filter(row -> row[:model] == "deterministic", tempDf2)[!, :cost]
                filteredPlaCost = filter(row -> row[:model] == inputString, tempDf2)[!, :cost]
                filteredRelativeCost = 100 * (filteredPlaCost./filteredDetCost)
                tempAverageRelativeCost[u, c] = mean(filteredRelativeCost)
                tempCostConfInterval[u, c] = 1.96 * std(filteredRelativeCost) / sqrt(length(filteredRelativeCost))
            end
        end
        # Store results in global matrix
        averageRelativeCost[:, (d-1)*2+1:d*2] = tempAverageRelativeCost
        costConfInterval[:, (d-1)*2+1:d*2] = tempCostConfInterval
    end

    return averageRelativeCost, costConfInterval
end
