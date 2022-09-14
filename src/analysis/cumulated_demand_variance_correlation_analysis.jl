function cumulated_demand_variance_correlation_analysis(forecastInput)
    ## cumulated_demand_variance_correlation_analysis(forecastInput)
    #
    # Input:
    #   - forecastInput
    # Output:
    #   - cumDemandVarianceAdd, cumDemandVarianceFw, cumDemandVarianceSamples

    forecast = forecastInput * ones(nbProducts, T)
    # Mean of forecast update
    meanForecastEvolAdd = zeros(nbProducts * T)
    meanVectorMult = [-multvar/2;-multvar/2]
    # Initialise matrices
    cumDemandVarianceAdd = zeros(nbCor)
    cumDemandVarianceFw = zeros(nbCor)
    cumDemandVarianceSamples = zeros(nbCor)

    # Calculate demand of cumulative demand in period 2 for each correlation coefficient
    for c in 1:nbCor
        ## Additive MMFE
        covarianceMatrix = [sigma1^2 sigma1 * sigma2 *correlationVector[c] ; sigma1 * sigma2 *correlationVector[c] sigma2^2]
        cumDisMatrix = cumulative_demand_distribution_matrix_from_additive_mmfe(forecast, covarianceMatrix, meanForecastEvolAdd , nbProducts)
        cumDemandVarianceAdd[c] = var(cumDisMatrix[t1, t2])
        ## Multiplicative MMFE
        covarianceMatrix = [multvar multvar * correlationVector[c] ; multvar *correlationVector[c] multvar]
        forecastDistribution = MvLogNormal(meanVectorMult, covarianceMatrix)
        # FW approximation
        cumDisMatrix = cumulative_demand_distribution_matrix_from_multiplicative_mmfe(forecast, meanVectorMult, covarianceMatrix, nbProducts)
        cumDemandVarianceFw[c] = var(cumDisMatrix[t1,t2])
        # Sampled
        correlatedSamples = rand(forecastDistribution, nbSamples)
        samplesDem1 = forecast[1,1] * vec(correlatedSamples[1, :])
        samplesDem2 = forecast[1,2] * (vec(correlatedSamples[2, :]) .* rand(LogNormal(-multvar/2, sqrt(multvar)), nbSamples))
        cumDemSamples = samplesDem1 + samplesDem2
        cumDemandVarianceSamples[c] = var(cumDemSamples)
    end

    return cumDemandVarianceAdd, cumDemandVarianceFw, cumDemandVarianceSamples
end
