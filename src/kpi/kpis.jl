function kpis(simInstanceList, backlogCost, demandRealStored, nbProducts, nbSimulations)
    ## kpis(simInstanceList, backlogCost, demandRealStored, nbProducts, nbSimulations)
    # Evaluate the performance of a simulation model in terms of cost and service level.
    # Input:
    #   - simInstanceList, backlogCost, demandRealStored, nbProducts, nbSimulations
    # Output:
    #   - simCosts: sum of inventory, setup and backlog costs
    #   - simSl: per-product service level
    #   - simAggSl: service level aggregated over all products

    simCosts = zeros(nbSimulations, 4)
    simSl = zeros(nbSimulations, nbProducts)
    simAggSl = zeros(nbSimulations)
    for s in 1:nbSimulations
        simCosts[s, :] = collect(realised_costs(backlogCost, simInstanceList[s]))
        simSl[s, :] = realised_service_level(demandRealStored[s], simInstanceList[s], nbProducts)[1]
        simAggSl[s] = realised_service_level(demandRealStored[s], simInstanceList[s], nbProducts)[2]
    end

    return simCosts, simSl, simAggSl
end

function kpis(simInstance :: SimulationInstance, backlogCost, demandRealStored, nbProducts)
    ## kpis(simInstance :: SimulationInstance, backlogCost, demandRealStored, nbProducts)
    # Evaluate the performance of a simulation model in terms of cost and service level.
    # Input:
    #   - simInstance :: SimulationInstance, backlogCost, demandRealStored, nbProducts
    # Output:
    #   - simCosts: sum of inventory, setup and backlog costs
    #   - simSl: per-product service level
    #   - simAggSl: service level aggregated over all products

    simCosts = collect(realised_costs(backlogCost, simInstance))
    simSl = realised_service_level(demandRealStored, simInstance, nbProducts)[1]
    simAggSl = realised_service_level(demandRealStored, simInstance, nbProducts)[2]

    return simCosts, simSl, simAggSl
end
