function extended_mmfe_with_pla(t, forecast, inventoryState, forecastEvolType, trueForecastEvolDis, capacity, tBreak, nbProducts; useFreeReturn=false)
    ## extended_mmfe_with_pla(t, forecast, inventoryState, forecastEvolType, trueForecastEvolDis, capacity, tBreak, nbProducts)
    # Build a scenario tree of demand from the forecast evolution, determine the cumulative demand distributions over the horizon, linearise the inventory and backlog functions, and solve the resulting PLA lot-sizing problem.
    # Input:
    #   - t, forecast, inventoryState, forecastEvolType, trueForecastEvolDis, capacity, tBreak, nbProducts
    # Output:
	#   - extendedPlaProdOut: Matrix[nbProducts, T], optimal production volumes
	#   - extendedProdRecourse: Matrix[nbProducts, T, nbScenarios], optimal recourse production volumes over scenario tree
	#   - extendedPlaSolveTime: scalar, solve time in seconds

	# --- Create multi-stage sceanrio tree from forecast evolution ---
    demandScenarios = build_scenario_tree_from_forecast_evolution(trueForecastEvolDis, forecast, nbProducts, forecastEvolType)

	# --- Obtain cumulative demand distribution ---
    if forecastEvolType == "Additive"
        cumulatedDemandDistributionMatrix = cumulative_demand_distribution_matrix_from_additive_mmfe(forecast, cov(trueForecastEvolDis), mean(trueForecastEvolDis), nbProducts)
    elseif forecastEvolType == "Multiplicative"
        cumulatedDemandDistributionMatrix =
        cumulative_demand_distribution_matrix_from_multiplicative_mmfe(forecast, location(trueForecastEvolDis), scale(trueForecastEvolDis), nbProducts)
    end
	# Linearise FOLF and c-FOLF over relevant domain
	uRange, folfOrigin, folfSlope, complFolfOrigin, complFolfSlope = folf_linearisation(inventoryState, cumulatedDemandDistributionMatrix, NB_SEGMENTS, capacity, nbProducts)

	# --- Solve optimisation problem ---
	if !(useFreeReturn)
	    extendedPlaProdOut, extendedProdRecourse, extendedPlaSolveTime = extended_pla_scenario_lot_sizing_with_prod_recourse(inventoryState, uRange, folfOrigin, folfSlope, complFolfOrigin, complFolfSlope, NB_SEGMENTS, backlogCost, capacity, nbProducts, demandScenarios, tBreak, gurobiEnv)
	elseif useFreeReturn
		# Linearise FOLF and c-FOLF over [inventoryState, 0] in first period to allow returns
		vRange, returnFolfOrigin, returnFolfSlope, returnComplFolfOrigin, returnComplFolfSlope = free_returns_folf_linearisation(inventoryState, cumulatedDemandDistributionMatrix, NB_RETURNS_SEGMENTS, capacity, nbProducts)
		check_return_linearization_consistency(vRange, returnFolfOrigin, folfOrigin, returnComplFolfOrigin, complFolfOrigin, nbProducts)
		# Solve lot-sizing model
	    extendedPlaProdOut, extendedProdRecourse, extendedPlaSolveTime = extended_pla_scenario_lot_sizing_free_returns_with_prod_recourse(inventoryState, uRange, folfOrigin, folfSlope, complFolfOrigin, complFolfSlope, vRange, returnFolfSlope, returnComplFolfSlope, NB_SEGMENTS, NB_RETURNS_SEGMENTS, backlogCost, capacity, nbProducts, demandScenarios, tBreak, gurobiEnv)
	end

    return extendedPlaProdOut, extendedProdRecourse, extendedPlaSolveTime
end
