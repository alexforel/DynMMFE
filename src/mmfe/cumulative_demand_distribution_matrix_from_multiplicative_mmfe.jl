function cumulative_demand_distribution_matrix_from_multiplicative_mmfe(forecast, logMeanUpdateVectorInput, logCovarMatrixInput, nbProducts)
    ## cumulative_demand_distribution_matrix_from_multiplicative_mmfe(forecast, logMeanUpdateVectorInput, logCovarMatrixInput, nbProducts)
    # Determine the cumulative demand distributions over the horizon from the input MMFE parameters using the Fenton-Wilkinson approximation.
    # Input:
    #   - forecast: Matrix[nbPRoducts, T], point-estimate forecast
    #   - logMeanUpdateVectorInput
    #   - logCovarMatrixInput
    #   - nbProducts
    # Output:
    #   - cumulatedDemandDistributionMatrix: Matrix[nbProducts, T ] of distributions

    # - Step 1 - Determine parameters of lognormal demand distribution in each period
    # 1.1/ Log-Mean vector of demand
    logMeanDemandVector = zeros(nbProducts * T)
    for k in 1:nbProducts
        for tau in 1:T
            indexToStore = (k-1) * T + tau
            ## Edge case: replace 0 or negative demand by very small demand
            demandForecast = forecast[k, tau]
            if demandForecast == 0
                demandForecast = 0.1
            end
            if demandForecast < 0
                @warn "demandForecast negative for prod $k at time $tau"
                demandForecast = 0.1
            end
            logMeanDemandVector[indexToStore] = log(demandForecast) + sum(logMeanUpdateVectorInput[i] for i in ((k-1) * T + 1):indexToStore)
        end
    end
    # 1.2/ Log-Covariance matrix of demand
    logCovarianceMatrixDemand = zeros(nbProducts * T, nbProducts * T)
    for k in 1:nbProducts
        for t1 in 1:T
            for t2 in 1:T
                logCovarianceMatrixDemand[t1 + (k-1) * T, t2 + (k-1) * T] = sum(
                    logCovarMatrixInput[t1-i+1 + (k-1) * T, t2-i+1 + (k-1) * T] for i in 1:min(t1,t2))
            end
        end
    end

    # - Step 2 - Approximate cumulative demand distribution using Fenton-Wilkinson approximation
    cumulatedDemandDistributionMatrix = lognormal_cumulative_demand_distribution_matrix(logMeanDemandVector, logCovarianceMatrixDemand, nbProducts, T)

    return cumulatedDemandDistributionMatrix
end

function lognormal_cumulative_demand_distribution_matrix(logMeanDemandVector, logCovarianceMatrixDemand, nbProducts, T)
    cumulatedDemandDistributionMatrix = Matrix(undef, nbProducts, T)
    for k in 1:nbProducts
        for tau in 1:T
            if tau == 1
                indexToRead = tau + (k-1) * T
                cumulatedDemandDistributionMatrix[k, tau] = LogNormal(logMeanDemandVector[indexToRead], sqrt(logCovarianceMatrixDemand[indexToRead, indexToRead]))
            else
                rangeOfDemandMatrices = (1 + (k-1) * T):(tau + (k-1) * T)
                cumulatedDemandDistributionMatrix[k, tau] = fenton_wilkinson_approx_sum_log_normal(logMeanDemandVector[rangeOfDemandMatrices], logCovarianceMatrixDemand[rangeOfDemandMatrices, rangeOfDemandMatrices])
            end
        end
    end

    return cumulatedDemandDistributionMatrix
end
