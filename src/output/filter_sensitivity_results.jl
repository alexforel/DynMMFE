function filter_sensitivity_results(inputDf)
    ## filter_sensitivity_results(inputDf)
    # Filter result data frame to the relevant sensitivity instance. Calculate confidence interval of relative costs and solution time. Perform hypothesis tests.
    # Input:
    #   - inputDf
    # Output:
    #   - averageRelativeCost, costConfInterval, solveTimeVector, timeConfInterval

    # Cost varying with tBreak
    averageRelativeCost = zeros(length(scenarioTreeVector), length(tBreakVector))
    costConfInterval = zeros(length(scenarioTreeVector), length(tBreakVector))
    for tree in 1:length(scenarioTreeVector)
        tempDf = filter(row -> row[:treeIndex] == tree, inputDf)
        for t in 1:length(tBreakVector)
            filteredRelativeCost = (filter(row -> row[:tbreak] == tBreakVector[t], tempDf)[!, :cost] ./ filter(row -> row[:tbreak] == 7, tempDf)[!, :cost]) .* 100
            averageRelativeCost[tree, t] = mean(filteredRelativeCost)
            costConfInterval[tree, t] = 1.96 * std(filteredRelativeCost) / sqrt(length(filteredRelativeCost))
        end
    end
    # Solve time varying with tBreak
    solveTimeVector = zeros(length(scenarioTreeVector), length(tBreakVector))
    timeConfInterval = zeros(length(scenarioTreeVector), length(tBreakVector))
    for tree in 1:length(scenarioTreeVector)
        tempDf = filter(row -> row[:treeIndex] == tree, inputDf)
        for t in 1:length(tBreakVector)
            filteredSolveTime = filter(row -> row[:tbreak] == tBreakVector[t], tempDf)[!, :solveTime]
            solveTimeVector[tree, t] = mean(filteredSolveTime)
            timeConfInterval[tree, t] = 1.96 * std(filteredSolveTime) / sqrt(length(filteredSolveTime))
        end
    end

    return averageRelativeCost, costConfInterval, solveTimeVector, timeConfInterval
end
