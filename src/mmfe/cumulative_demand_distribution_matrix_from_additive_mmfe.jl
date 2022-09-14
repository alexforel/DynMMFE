function cumulative_demand_distribution_matrix_from_additive_mmfe(forecast, covarianceMatrixForecastEvolVector, meanForecastEvolVector, nbProducts)
    ## cumulative_demand_distribution_matrix_from_additive_mmfe(forecast, covarianceMatrixForecastEvolVector, meanForecastEvolVector, nbProducts)
    # Determine the cumulative demand distributions over the horizon from the input MMFE parameters.
    # Input:
    #   - forecast: Matrix[nbPRoducts, T], point-estimate forecast
    #   - covarianceMatrixForecastEvolVector
    #   - meanForecastEvolVector
    #   - nbProducts
    # Output:
    #   - cumulatedDemandDistributionMatrix: Matrix[nbProducts, T ] of distributions

    # Determine the parameters of the cumulative demand distribution for each product and time period
    meanDemandVector = zeros(nbProducts, T)
    demandCovarianceMatrix = zeros(nbProducts, T, T)
    for k in 1:nbProducts
        # Step 1 - Determine demand distribution from forecast evolution
        # Mean
        for t in 1:T
            # Mean of demand in each period = forecast + update biases
            meanDemandVector[k, t] = forecast[k, t] + sum(meanForecastEvolVector[((k-1)*T+1):((k-1)*T+t)])
        end
        # Variance
        productCovarMatrix = covarianceMatrixForecastEvolVector[((k-1)*T+1):((k-1)*T+T), ((k-1)*T+1):((k-1)*T+T)]
        for t1 in 1:T
            for t2 in 1:T
                demandCovarianceMatrix[k, t1, t2] = sum(productCovarMatrix[t1 - tau + 1, t2 - tau + 1] for tau in 1:min(t1, t2))
            end
        end
    end

    # Step 2 - Determine cumulative demand distribution
    cumulatedDemandDistributionMatrix = normal_cumulative_demand_distribution_matrix(meanDemandVector, demandCovarianceMatrix, nbProducts, T)

    return cumulatedDemandDistributionMatrix
end

function normal_cumulative_demand_distribution_matrix(meanDemandVector, demandCovarianceMatrix, nbProducts, T)
    # Initialise Matrix of distributions
    cumulatedDemandDistributionMatrix = Matrix(undef, nbProducts, T)
    for k in 1:nbProducts
        meanCumulatedDemand = zeros(T)
        varCumulativeDemand = zeros(T)
        for t in 1:T
            # Mean of cumulative demand
            meanCumulatedDemand[t] = sum(meanDemandVector[k, 1:t])
            # Variance of cumulative demand
            varCumulativeDemand[t] = sum(sum(demandCovarianceMatrix[k, t1, t2] for t1 in 1:t) for t2 in 1:t)
            # Create distribution and store it in matrix
            if (meanCumulatedDemand[t] - sqrt(varCumulativeDemand[t]))  < 0
                # Create a truncated distribution if probability of negative demand is too high
                cumulatedDemandDistributionMatrix[k, t] = truncated(Normal(meanCumulatedDemand[t], sqrt(varCumulativeDemand[t])), 0.0, Inf)
            else
                cumulatedDemandDistributionMatrix[k, t] = Normal(meanCumulatedDemand[t], sqrt(varCumulativeDemand[t]))
            end
        end
    end
    return cumulatedDemandDistributionMatrix
end
