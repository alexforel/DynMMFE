function demand_driven_stochastic(t, forecast, stdMatrix, inventoryState, capacity, nbProducts, forecastEvolType)
    ## demand_driven_stochastic(t, forecast, stdMatrix, inventoryState, capacity, nbProducts, forecastEvolType)
    # Determine cumulative demand distributions over horizon from demand distributions and solve stochastic lot sizing model with PLA.
    # Input:
    #   - t, forecast, stdMatrix, inventoryState, capacity, nbProducts, forecastEvolType
    # Output:
    #   - plaProdOut, plaSolveTime

    # Cumulative demand distribution matrix
    if forecastEvolType == "Additive"
        meanDemVec = forecast
        demandCovarianceMatrix = zeros(nbProducts, T, T)
        for k in 1:nbProducts
            demandCovarianceMatrix[k, :, :] = diagm(stdMatrix[k, :].^2)
        end
        cumulatedDemandDistributionMatrix = normal_cumulative_demand_distribution_matrix(meanDemVec, demandCovarianceMatrix, nbProducts, T)
    elseif forecastEvolType == "Multiplicative"
        logMeanDemandVec = vec(forecast')
        logCovarianceMatrixDemand = diagm(vec((stdMatrix.^2)'))
        cumulatedDemandDistributionMatrix = lognormal_cumulative_demand_distribution_matrix(logMeanDemandVec, logCovarianceMatrixDemand, nbProducts, T)
    end

    # Linearise FOLF and c-FOLF over relevant domain
    uRange, folfOrigin, folfSlope, complFolfOrigin, complFolfSlope = folf_linearisation(inventoryState, cumulatedDemandDistributionMatrix, NB_SEGMENTS, capacity, nbProducts)

    # --- Solve optimisation problem ---
    plaProdOut, plaSolveTime = pla_lot_sizing_no_prod_recourse(inventoryState, uRange, folfOrigin, folfSlope, complFolfOrigin, complFolfSlope, NB_SEGMENTS, backlogCost, capacity, nbProducts, gurobiEnv)

    return plaProdOut, plaSolveTime
end
