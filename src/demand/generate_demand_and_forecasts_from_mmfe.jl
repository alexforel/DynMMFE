function generate_demand_and_forecasts_from_mmfe(trueForecastEvolDis, demandDistributions, T, NB_SIMULATIONS, forecastEvolType)
    ## generate_demand_and_forecasts_from_mmfe(trueForecastEvolDis, demandDistributions, T, NB_SIMULATIONS, forecastEvolType)
    # Simulate forecast evolution for each sim over the simulation length. The forecast evolution also describes the demand realisation. A new forecast is added to the horizon following the demand pattern.
    # Input:
    #   - trueForecastEvolDis: MvDistribution, multivariate forecast update distribution
    #   - demandDistributions: Array of Distributions, different mean and std for each period
    #   - T, NB_SIMULATIONS, forecastEvolType
    # Output:
    #   - demandRealStored, forecastStored

    # Initialise demand related arrays to store data
    demandRealStored, forecastStored = init_simulation_demand_matrix(NB_SIMULATIONS, simulationLength, NB_PRODUCTS, T)

    # Create demand realisation and forecasts over simulation length
    for simIndex in 1:NB_SIMULATIONS
        # First demand forecast
        for t in 1:T
            forecastStored[simIndex][:, t, 1] = max.(rand(demandDistributions[t], NB_PRODUCTS), 0.1) # forbid negative and zero demand/forecast
        end
        # Remove demand in first periods if warm-up phase
        if preSimT > 0
            forecastStored[simIndex][:, 1:preSimT, 1] = zeros(NB_PRODUCTS, preSimT)
        end

        for t in 1:simulationLength
            # Update forecast and demand realisation
            update_demand_and_forecast!(t, trueForecastEvolDis, forecastStored[simIndex], demandRealStored[simIndex], forecastEvolType, NB_PRODUCTS)
            # Add new forecast
            lastPeriodTimingInSimulation = t + 1 + T - 1
            if lastPeriodTimingInSimulation <= (preSimT + inSimT )
                forecastStored[simIndex][:, T, t+1] = max.(rand(demandDistributions[lastPeriodTimingInSimulation], NB_PRODUCTS), 0.1) # forbid negative and zero demand/forecast
            end
        end
    end

    # Output
    return demandRealStored, forecastStored
end

function init_simulation_demand_matrix(NB_SIMULATIONS, simulationLength, NB_PRODUCTS, T)
    ## init_simulation_demand_matrix(NB_SIMULATIONS, simulationLength, NB_PRODUCTS, T)
    # Initialise the demand and forecasts arrays with matrices of zeros.
    # Input:
    #   - NB_SIMULATIONS, simulationLength, NB_PRODUCTS, T
    # Output:
    #   - demandRealStored, forecastStored

    # Create arrays to store simulationInstances
    demandRealStored = Matrix(undef, NB_SIMULATIONS, 1)
    forecastStored = Matrix(undef, NB_SIMULATIONS, 1)
    for simIndex in 1:NB_SIMULATIONS
        demandRealStored[simIndex] = zeros(NB_PRODUCTS, simulationLength)
        forecastStored[simIndex] = zeros(NB_PRODUCTS, T, simulationLength)
    end

    return demandRealStored, forecastStored
end
