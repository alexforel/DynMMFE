function realised_costs(backlogCost, simInstance :: SimulationInstance)
    ## realised_costs(backlogCost, simInstance :: SimulationInstance)
    # Calculate the inventory, setup and backordering costs and their sum over the simulation length.
    # Input:
    #   - backlogCost, simInstance :: SimulationInstance
    # Output:
    #   - invCostSimInstance, backlogCostSimInstance, setupCostSimInstance, totalCost

    # Inventory costs
    invCostSimInstance = sum(invCost' * max.(simInstance.inventoryState, 0))
    # Backordering costs
    backlogCostSimInstance = sum(backlogCost' * max.(-simInstance.inventoryState, 0))
    # Setup costs
    setupCostSimInstance = sum(setupCost * simInstance.isSetup)
    # Total costs = sum of all costs
    totalCost = invCostSimInstance + backlogCostSimInstance + setupCostSimInstance

    return invCostSimInstance, backlogCostSimInstance, setupCostSimInstance, totalCost
end
