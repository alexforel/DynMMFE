function true_forecast_evolution_and_estimated_distribution(varUpdate, demandDistributions, timeCorrelation, productCorrelation, forecastEvolType ; useDifferentProductVars = false)
    ## true_forecast_evolution_and_estimated_distribution(varUpdate, demandDistributions, timeCorrelation, productCorrelation, forecastEvolType)
    # Define the multivariate forecast evolution distribution and simulate the corresponding "mismatched" distribution.
    # Input:
    #   - varUpdate, demandDistributions, timeCorrelation, productCorrelation, forecastEvolType
    # Output:
    #   - multForecastEvolDis, multMeanVector, addMmfeMean, addMmfeCovarianceMatrix, ddsMean, ddsStd

    NB_SIM = Int(1e6)

    # Define the true multivariate forecast evolution process from input parameters
    if forecastEvolType == "Additive"
        addForecastEvolDis = create_forecast_multivariate_distributions(varUpdate, timeCorrelation, productCorrelation, NB_PRODUCTS, T, forecastEvolType, useDifferentProductVars)
    elseif forecastEvolType == "Multiplicative"
        # True model
        multForecastEvolDis = create_forecast_multivariate_distributions(varUpdate, timeCorrelation, productCorrelation, NB_PRODUCTS, T, forecastEvolType, useDifferentProductVars)
    else
        printn("Forecast evolution model specified is incorrect: $forecastEvolType.")
        printn("Should be either \"Additive\" or \"Multiplicative\".")
    end

    ## ---- Simulate the forecast evolution process and measure covariance matrix of mismatched MMFE update ----
    forecastHistory = zeros(NB_PRODUCTS * T, NB_SIM + 1)
    postUpdateForecast = zeros(NB_PRODUCTS * T, NB_SIM)
    mismatchForecastUpdateHistory = zeros(NB_PRODUCTS * T, NB_SIM)
    # Initial forecast in first period
    for k in 1:NB_PRODUCTS
        for t in 1:T
            forecastHistory[T * (k-1) + t, 1] = rand(demandDistributions[t])
        end
    end
    # Simulate over NB_SIM periods
    for s in 1:NB_SIM
        # Sample demand realisation and forecast update
        if forecastEvolType == "Additive"
            sampleUpdate = rand(addForecastEvolDis, 1)
        elseif forecastEvolType == "Multiplicative"
            sampleUpdate = rand(multForecastEvolDis, 1)
        end
        # Store new forecast and measure the mismatched forecast update
        if forecastEvolType == "Additive"
            postUpdateForecast[:, s] = max.(forecastHistory[:, s] + sampleUpdate, 0.0001) # Forbid negative forecasts
            mismatchForecastUpdateHistory[:, s] = log.(postUpdateForecast[:, s] ./ forecastHistory[:, s])
        elseif forecastEvolType == "Multiplicative"
            postUpdateForecast[:, s] = forecastHistory[:, s] .* sampleUpdate
            mismatchForecastUpdateHistory[:, s] = postUpdateForecast[:, s] - forecastHistory[:, s]
        end
        # Update forecast over common horizon and add new
        for k in 1:NB_PRODUCTS
            sInd = T * (k-1) + 1
            eInd =  T * k - 1
            if forecastEvolType == "Additive"
                forecastHistory[sInd:eInd, s+1] = max.(sampleUpdate[(sInd+1):(eInd+1)] + forecastHistory[(sInd+1):(eInd+1), s], 0.0001) # Forbid negative forecast
            elseif forecastEvolType == "Multiplicative"
                forecastHistory[sInd:eInd, s+1] = sampleUpdate[(sInd+1):(eInd+1)] .* forecastHistory[(sInd+1):(eInd+1), s]
            end
            # New forecast for the last period in horizon
            if ((s+1) % simulationLength) == 0
                periodToSample = simulationLength
            else
                periodToSample = ((s+1) % simulationLength)
            end
            forecastHistory[T * k, s+1] = rand(demandDistributions[periodToSample])
        end
    end

    ## Output
    if forecastEvolType == "Additive"
        addMmfeMean = mean(addForecastEvolDis)
        addMmfeCovarianceMatrix = cov(addForecastEvolDis)
        # Gather indices of all simulations to remove because of 0 forecast value
        indicesToRemove = []
        for s in 1:NB_SIM
            # Remove the forecast vector if any of its element is smaller than the tolerance
            if sum(postUpdateForecast[:, s] .<= 1) >= 1
                # Remove the forecast update at index s because it leads to a null value
                append!(indicesToRemove, s)
                # Remove the forecast update at index s+1 because it updates a null value
                append!(indicesToRemove, s+1)
            end
        end
        indicesForEstimation = setdiff(1:NB_SIM, unique!(indicesToRemove))
        percentageRemoved = round(100 - length(indicesForEstimation)/NB_SIM * 100)
        println("$forecastEvolType , var: $varUpdate : Removed $percentageRemoved percent of update vectors.")
        ## Determine distribution
        multMeanVector = mean(mismatchForecastUpdateHistory[:, indicesForEstimation], dims = 2)
        multCovarianceMatrix = cov(mismatchForecastUpdateHistory[:, indicesForEstimation]')
        multForecastEvolDis = MvLogNormal(vec(multMeanVector), multCovarianceMatrix)
    elseif forecastEvolType == "Multiplicative"
        addMmfeMean = mean(mismatchForecastUpdateHistory, dims = 2)
        addMmfeCovarianceMatrix = cov(mismatchForecastUpdateHistory')
        multMeanVector = location(multForecastEvolDis)
    end

    ## ----- Demand-driven stochastic benchmark -----
    # Take into account demand seasonality
    # Mean = average of demand of each product / period
    # Std = standard deviation of demand of each product / period
    ddsMean = zeros(NB_PRODUCTS, simulationLength)
    ddsStd = zeros(NB_PRODUCTS, simulationLength)
    for k in 1:NB_PRODUCTS
        for i in 1:simulationLength
            if (i-1) == simulationLength
                indices = (collect(1:NB_SIM) .% simulationLength) .== 0
            else
                indices = (collect(1:NB_SIM) .% simulationLength) .== (i-1)
            end
            if forecastEvolType == "Additive"
                relevantSamples = postUpdateForecast[1 + (k-1) * T, indices]
            elseif forecastEvolType == "Multiplicative"
                relevantSamples = log.(postUpdateForecast[1 + (k-1) * T, indices])
            end
            ddsMean[k, i] = mean(relevantSamples)
            ddsStd[k, i] = std(relevantSamples)
        end
    end

    return multForecastEvolDis, multMeanVector, addMmfeMean, addMmfeCovarianceMatrix, ddsMean, ddsStd
end

function create_forecast_multivariate_distributions(varUpdate, timeCorrelation, productCorrelation, NB_PRODUCTS, T, forecastEvolType, useDifferentProductVars)
    ## create_forecast_multivariate_distributions(timeCorrelation, NB_PRODUCTS, T; isCorrelated = true)
    # Returns the multivariate forecast evolution distribution
    # Input:
    #   - varUpdate, timeCorrelation, productCorrelation, NB_PRODUCTS, T
    #   - forecastEvolType: String, either Additive or Multiplicative MMFE
    # Output:
    #   - forecastEvolDistribution: MultivariateDistribution, of dimension NB_PRODUCTS times T

    # Create covariance matrix from variance vector
    singleProductTimeCovariance = diagm(varUpdate)
    if useDifferentProductVars
        # Products are uncorrelated
        prod2VarCoeff = 7 / 4
        updateCovarMatrix = [singleProductTimeCovariance singleProductTimeCovariance * 0 ; singleProductTimeCovariance * 0 prod2VarCoeff*singleProductTimeCovariance]
    else
        # We only consider correlation between time periods 1 and 2
        singleProductTimeCovariance[1, 2] = varUpdate[1] * timeCorrelation
        singleProductTimeCovariance[2, 1] = varUpdate[1] * timeCorrelation
        productCovarMatrix = singleProductTimeCovariance * productCorrelation
        updateCovarMatrix = [singleProductTimeCovariance productCovarMatrix ; productCovarMatrix singleProductTimeCovariance] # ! only valid for two products
    end

    # Create multivariate distribution objects
    if forecastEvolType == "Additive"
        forecastEvolDistribution = MvNormal(zeros(NB_PRODUCTS * T), updateCovarMatrix)
    elseif forecastEvolType == "Multiplicative"
        forecastEvolDistribution = MvLogNormal(-0.5 * diag(updateCovarMatrix), updateCovarMatrix)
    end

    ## Output
    return forecastEvolDistribution
end
