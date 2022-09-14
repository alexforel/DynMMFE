struct SimulationInstance
  ## SimulationInstance(bomDf, productIdList)

  # Decisions
  referencePlan
  productionPlanTree   # Only for stochastic models
  inventoryPlan
  backorderPlan
  # State
  inventoryState
  productionRealisation
  isSetup
  # Solver related
  solveTime

  function SimulationInstance(initInv, simulationLength, nbProducts)
    # Input:
    #   - initInv
    #   - simulationLength
    #   - nbProducts
    # Output:
    #   - productNameList

    referencePlan = Matrix(undef, simulationLength, 1)
    productionPlanTree = Matrix(undef, simulationLength, 1)
    inventoryPlan  = Matrix(undef, simulationLength, 1)
    backorderPlan  = Matrix(undef, simulationLength, 1)
    # State
    inventoryState  = zeros(nbProducts, simulationLength + 1)
    productionRealisation = zeros(nbProducts, simulationLength)
    isSetup = zeros(nbProducts, simulationLength)
    solveTime = zeros(simulationLength)
    # Initial conditions
    inventoryState[:, 1] = initInv

    new(referencePlan, productionPlanTree, inventoryPlan, backorderPlan, inventoryState, productionRealisation, isSetup, solveTime)
  end
end
