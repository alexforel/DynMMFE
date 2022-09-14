function update_instance_state!(t, prodPlan, solveTime, demandRealisation, simInstance, nbProducts)
    ## update_instance_state!(t, prodPlan, solveTime, demandRealisation, simInstance, nbProducts)
    # Update the state stored in the simulationInstance object.
    # Input:
    #   - t, prodPlan, solveTime, demandRealisation, simInstance, nbProducts
    # Output:
    #   - nothing

    simInstance.referencePlan[t] = prodPlan
    simInstance.productionRealisation[:, t] = prodPlan[:, 1]
    # Update state using inventory balance equation: add production and substract demand
    simInstance.inventoryState[:, t + 1] = simInstance.inventoryState[:, t] + simInstance.productionRealisation[:, t] - demandRealisation
    for k in 1:nbProducts
        if simInstance.productionRealisation[k, t] > 1e-4
                simInstance.isSetup[k, t] = 1
        end
    end
    simInstance.solveTime[t] = solveTime

    return nothing
end

function update_instance_state!(t, prodPlan, prodRecourse, solveTime, demandRealisation, simInstance, nbProducts)
    ## update_instance_state!(t, prodPlan, prodRecourse, solveTime, demandRealisation, simInstance, nbProducts)
    # Function for stochastic model with recourse

    update_instance_state!(t, prodPlan, solveTime, demandRealisation, simInstance, nbProducts)
    simInstance.productionPlanTree[t] = prodRecourse
    return nothing
end
