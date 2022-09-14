function build_scenario_tree_from_forecast_evolution(updateDistribution, forecast, nbProducts, modelTypeString; monteCarloSampling = false)
    ## build_scenario_tree_from_forecast_evolution(updateDistribution, forecast, nbProducts, modelTypeString; monteCarloSampling = false)
    # Create multi-stage scenario tree from forecast evolution distribution. Multidimension latin hypercube strate are used to sample the multivariate forecast evolution process
    # Input:
    #   - updateDistribution, forecast, nbProducts, modelTypeString; monteCarloSampling = false
    # Output:
    #   - demandScenarios: Matrix [nbProducts, T, NB_SCENARIOS], demand scenarios

    # Initialise scenario matrix
    demandScenarios = zeros(nbProducts, T, NB_SCENARIOS)
    for s in 1:NB_SCENARIOS
        demandScenarios[:, :, s] = forecast
    end
    # Generate updates from period 1 to T and add/multiply them to forecasts following MMFE
    for t in 1:T
        nbNodes = NB_SCEN_NODES[t]
        # Sample from the latin hypercube strata
        if modelTypeString == "Additive"
            latinSamples = latin_hypercube_sampling_with_multidimensional_uniformity(LATIN_STRATA[t], updateDistribution, nbNodes)
        elseif modelTypeString == "Multiplicative"
            tLoc = location(updateDistribution)
            tScale = scale(updateDistribution)
            normDis = MvNormal(tLoc, tScale)
            latinSamples = exp.(latin_hypercube_sampling_with_multidimensional_uniformity(LATIN_STRATA[t], normDis, nbNodes))
        end
        # Update the initial forecast following the MMFE process over the prediction horizon
        for node in 1:nbNodes
            if monteCarloSampling
                update = rand(updateDistribution)
            else
                update = latinSamples[:, node]
            end
            lengthBranch = Int(NB_SCENARIOS / nbNodes)
            for s in ((node-1) * lengthBranch + 1):(node * lengthBranch)
                for k in 1:nbProducts
                    for tau in 1:(T-t+1)
                        demandTimeIndex = tau + t - 1
                        if modelTypeString == "Additive"
                            demandScenarios[k, demandTimeIndex, s] = demandScenarios[k, demandTimeIndex, s] + update[(k-1)*T + tau]
                        elseif modelTypeString == "Multiplicative"
                            demandScenarios[k, demandTimeIndex, s] = demandScenarios[k, demandTimeIndex, s] * update[(k-1)*T + tau]
                        end
                    end
                end
            end
        end
    end
    demandScenarios = max.(demandScenarios, 0.0)

    return demandScenarios
end
