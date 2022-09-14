function update_demand_and_forecast!(t, forecastEvolutionDistribution, forecastStoredInstance, demandRealStoredInstance, forecastEvolType, nbProducts)
    ## update_demand_and_forecast!(forEvolDisArray, forecastStored, demandRealStored)
    # Sample a forecast updater from MMFE distribution and updates the forecasts/demand with it.
    # Input:
    #   - t, forecastEvolutionDistribution, forecastStoredInstance, demandRealStoredInstance, forecastEvolType, nbProducts
    # Output:
    #   nothing

    # Sample demand realisation and forecast update from MMFE distribution
    update = rand(forecastEvolutionDistribution, 1)

    # -- Forecast update --
    # Reformat forecast update in matrix form
    updateReshape = zeros(nbProducts, T)
    for k in 1:nbProducts
        updateReshape[k, :] = update[((k-1) * T + 1):(k * T), 1]
    end
    # Update forecast to demand
    if forecastEvolType == "Additive"
        demandRealStoredInstance[:, t] = forecastStoredInstance[:, 1, t] + updateReshape[:, 1]
    elseif forecastEvolType == "Multiplicative"
        demandRealStoredInstance[:, t] = forecastStoredInstance[:, 1, t] .* updateReshape[:, 1]
    end
    # Only positive values allowed
    demandRealStoredInstance[:, t] = max.(demandRealStoredInstance[:, t], 0)
    # Update forecast vector ove rhorizon
    if t < simulationLength
        if forecastEvolType == "Additive"
            forecastStoredInstance[:, 1:(T-1), t+1] = forecastStoredInstance[:, 2:T, t] + updateReshape[:, 2:T]
        elseif forecastEvolType == "Multiplicative"
            forecastStoredInstance[:, 1:(T-1), t+1] = forecastStoredInstance[:, 2:T, t] .* updateReshape[:, 2:T]
        end
        # Only positive values allowed
        forecastStoredInstance[:, 1:(T-1), t+1] = max.(forecastStoredInstance[:, 1:(T-1), t+1], 0)

        # No demand after end of simulation
        for tau in 1:T
            if ((t + 1 + tau - 1) > (preSimT + inSimT))
                forecastStoredInstance[:, tau, t+1] .= 0
            end
        end
    end

    return nothing
end
