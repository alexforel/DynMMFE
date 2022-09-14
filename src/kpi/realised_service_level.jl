function realised_service_level(demandRealisation, simInstance, nbProducts)
    ## realised_service_level(demandRealisation, simInstance, nbProducts)
    # Evaluate the acheived service level per product and aggregated.
    # Input:
    #   - demandRealisation, simInstance, nbProducts
    # Output:
    #   - realisedServiceLevel, aggSl

    realisedServiceLevel = zeros(nbProducts)
    # Achieved service level per product
    for p in 1:nbProducts
        sumDemandHorizon = sum(demandRealisation[p, :])
        if sumDemandHorizon == 0
            realisedServiceLevel[p] = 1
        else
            realisedServiceLevel[p] = 1 - (sum(max.(-simInstance.inventoryState[p, :], 0)) / sumDemandHorizon)
        end
    end
    # Aggregated service level over all products
    aggSl = sum(realisedServiceLevel[p] * sum(demandRealisation[p, :]) / sum(demandRealisation) for p in 1:nbProducts)

    return realisedServiceLevel, aggSl
end
