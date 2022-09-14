function additive_mmfe_with_pla(t, forecast, inventoryState, capacity, nbProds, addMeanVector, addMmfeCovarianceMatrix; useFreeReturn=false)
    ## additive_mmfe_with_pla(t, forecast, inventoryState, capacity, nbProds, addMeanVector, addMmfeCovarianceMatrix)
    # Determine the cumulative demand distributions over the horizon, linearise the inventory and backlog functions, and solve the resulting PLA lot-sizing problem.
    # Input:
    #   - t, forecast, inventoryState, capacity, nbProds, addMeanVector, addMmfeCovarianceMatrix
    # Output:
	#   - plaProdOut: Matrix[nbProducts, T], optimal production volumes
	#   - plaSolveTime: scalar, solve time in seconds

    # Determine distribution of each cumulative demand
    cumulatedDemandDistributionMatrix = cumulative_demand_distribution_matrix_from_additive_mmfe(forecast, addMmfeCovarianceMatrix, addMeanVector, nbProds)

    # Linearise FOLF and c-FOLF over relevant domain
    uRange, folfOrigin, folfSlope, complFolfOrigin, complFolfSlope = folf_linearisation(inventoryState, cumulatedDemandDistributionMatrix, NB_SEGMENTS, capacity, nbProds)

	# --- Solve optimisation problem ---
	if !(useFreeReturn)
    	plaProdOut, plaSolveTime = pla_lot_sizing_no_prod_recourse(inventoryState, uRange, folfOrigin, folfSlope, complFolfOrigin, complFolfSlope, NB_SEGMENTS, backlogCost, capacity, nbProds, gurobiEnv)
	elseif useFreeReturn
		# Linearise FOLF and c-FOLF over [inventoryState, 0] in first period to allow returns
		vRange, returnFolfOrigin, returnFolfSlope, returnComplFolfOrigin, returnComplFolfSlope = free_returns_folf_linearisation(inventoryState, cumulatedDemandDistributionMatrix, NB_RETURNS_SEGMENTS, capacity, nbProds)
		check_return_linearization_consistency(vRange, returnFolfOrigin, folfOrigin, returnComplFolfOrigin, complFolfOrigin, nbProds)
		# Solve lot-sizing model
	    plaProdOut, plaSolveTime = pla_lot_sizing_free_returns_no_prod_recourse(inventoryState, uRange, folfOrigin, folfSlope, complFolfOrigin, complFolfSlope, vRange, returnFolfSlope, returnComplFolfSlope, NB_SEGMENTS, NB_RETURNS_SEGMENTS, backlogCost, capacity, nbProds, gurobiEnv)
	end

    return plaProdOut, plaSolveTime
end
