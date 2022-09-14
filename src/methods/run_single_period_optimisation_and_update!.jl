function run_single_period_optimisation_and_update!(t, forecast, demand, detSimInstance, addMmfeSimInstance, addMmfeMean, addMmfeCovarianceMatrix, multMmfeSimInstance, multForecastEvolDis, multMeanVector, extendMmfeSimInstance, forecastEvolType, trueForecastEvolDis, capacity, tBreak, nbProducts)
    ## run_single_period_optimisation_and_update!(t, forecast, demand, detSimInstance, addMmfeSimInstance, addMmfeMean, addMmfeCovarianceMatrix, multMmfeSimInstance, multForecastEvolDis, multMeanVector, extendMmfeSimInstance, forecastEvolType, trueForecastEvolDis, capacity, tBreak, nbProducts)
    # Solve optimisation problems of all models and update state with production decision and observed demand.
    # Input:
    #   - t, forecast, demand, detSimInstance, addMmfeSimInstance, addMmfeMean, addMmfeCovarianceMatrix, multMmfeSimInstance, multForecastEvolDis, multMeanVector, extendMmfeSimInstance, forecastEvolType, trueForecastEvolDis, capacity, tBreak, nbProducts
    # Output:
    #   - nothing

    # Solve deterministic model
    detProdOut, detSolveTime = deterministic_equivalent_optim(detSimInstance.inventoryState[:, t], forecast, capacity, nbProducts, gurobiEnv)
    # Solve problem using PLA
    addProdOut, addPlaSolveTime = additive_mmfe_with_pla(t, forecast, addMmfeSimInstance.inventoryState[:, t], capacity, nbProducts, addMmfeMean, addMmfeCovarianceMatrix)
    multProdOut, multPlaSolveTime = multiplicative_mmfe_with_pla(t, forecast, multMmfeSimInstance.inventoryState[:, t], multMeanVector, scale(multForecastEvolDis), capacity, nbProducts)
    extendPlaProdOut, extendProdRecourse, extendPlaSolveTime = extended_mmfe_with_pla(t, forecast, extendMmfeSimInstance.inventoryState[:, t], forecastEvolType, trueForecastEvolDis, capacity, tBreak, nbProducts)
    # Update simulation instances
    update_instance_state!(t, detProdOut, detSolveTime, demand, detSimInstance, nbProducts)
    update_instance_state!(t, addProdOut, addPlaSolveTime, demand, addMmfeSimInstance, nbProducts)
    update_instance_state!(t, multProdOut, multPlaSolveTime, demand, multMmfeSimInstance, nbProducts)
    update_instance_state!(t, extendPlaProdOut, extendProdRecourse, extendPlaSolveTime, demand, extendMmfeSimInstance, nbProducts)

    # Output
    return nothing
end
