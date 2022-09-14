function init_simulation_instances(NB_SIMULATIONS, INIT_INV, simulationLength, NB_PRODUCTS, T)
    ## init_simulation_instances(NB_SIMULATIONS, INIT_INV, simulationLength, NB_PRODUCTS, T)
    # Initialise the matrix of simulation instances
    # Input:
    #   - NB_SIMULATIONS, INIT_INV, simulationLength, NB_PRODUCTS, T
    # Output:
    #   - simInstanceMatrix Matrix[NB_SIMULATIONS], each element in a simulation instance

    # Create arrays to store simulationInstances
    simInstanceMatrix = Matrix(undef, NB_SIMULATIONS, 1)
    for simIndex in 1:NB_SIMULATIONS
        simInstanceMatrix[simIndex] = SimulationInstance(INIT_INV, simulationLength, NB_PRODUCTS)
    end


    return simInstanceMatrix
end
