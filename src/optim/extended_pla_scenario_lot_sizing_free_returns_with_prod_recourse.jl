function extended_pla_scenario_lot_sizing_free_returns_with_prod_recourse(initInv, uRange, folfOrigin, folfSlope, complFolfOrigin, complFolfSlope, vRange, returnFolfSlope, returnComplFolfSlope, nbSegments, nbReturnsSegments, backlogCost, capacity, nbProducts, demandScenarios, tBreak, gurobiEnv)
	## extended_pla_scenario_lot_sizing_free_returns_with_prod_recourse(initInv, uRange, folfOrigin, folfSlope, complFolfOrigin, complFolfSlope, nbSegments, backlogCost, capacity, nbProducts, demandScenarios, tBreak, gurobiEnv)
	# Formulate a stochastic lot-sizing with piecewise-linearisation of the inventory and backlog functions and a scenario-based recourse formulation, and solve it using Gurobi.
	# Input:
	#   - initInv, uRange, folfOrigin, folfSlope, complFolfOrigin, complFolfSlope, nbSegments, backlogCost, capacity, nbProducts, demandScenarios, tBreak, gurobiEnv
	# Output:
	#   - prodOut, productionRecourseOut, solveTime

	if VERBOSE == 1
		println("---- Solving Extended PLA-SCEN model ----")
	end
	## Define Optimisation model
	lotSizingModel = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(gurobiEnv), "TimeLimit" => 100, "MIPGap" => 0.01, "OutputFlag" => VERBOSE))
	# Read parameters from inputs
	nbProductionNodes = NB_SCEN_NODES[(tBreak-1):(T-1)]
	nbInventoryNodes = NB_SCEN_NODES
	nbScenarios = nbInventoryNodes[end]
	## ---- Variables ----
	@variables(lotSizingModel,
		begin
			# Production variables
			production[1:nbProducts, 1:(tBreak - 1)] >= 0
			productionRecourse[1:nbProducts, 1:(T - tBreak + 1), 1:nbScenarios] >= 0
			setup[1:nbProducts, 1:T], Bin
			# Inventory variables
			expInv[1:nbProducts, 1:(tBreak - 1)]
			expBacklog[1:nbProducts, 1:(tBreak - 1)]
			invRecourse[1:nbProducts, 1:T, 1:nbScenarios] >= 0
			backorderRecourse[1:nbProducts, 1:T, 1:nbScenarios] >= 0
			# Auxiliary variables
			lambda[1:nbProducts, 1:(tBreak - 1), 1:nbSegments], Bin
			wSegment[1:nbProducts, 1:(tBreak - 1), 1:nbSegments] >= 0
			## Free return variables
			returnVol[1:nbProducts] >= 0
			# Auxiliary variables
			nu[1:nbProducts, 1:nbReturnsSegments], Bin
			omegaSegment[1:nbProducts, 1:nbReturnsSegments] >= 0
		end)

	## ---- Objective ----
	# Objective: min production cost, inventory costs and shortfall costs
	plaCost = sum(sum(invCost[k] * expInv[k,t] + backlogCost[k] * expBacklog[k,t]  for t in 1:(tBreak-1)) for k in 1:nbProducts)
	scenCost = 1/nbScenarios * sum(sum(invCost' * invRecourse[:, t, s] + backlogCost' * backorderRecourse[:, t, s] for t in tBreak:T) for s in 1:nbScenarios)
	setCost = setupCost * sum(setup)
	@objective(lotSizingModel, Min, plaCost + scenCost + setCost)

	## ---- Constraints ----
	# Calculate expected inventory and backlog
	@constraint(lotSizingModel, [k in 1:nbProducts, t in 1:(tBreak - 1)],
		expInv[k, t] == complFolfOrigin[k, t] + sum(complFolfSlope[k, t, l] * wSegment[k, t, l] for l in 1:nbSegments) - sum(returnComplFolfSlope[k, l] * omegaSegment[k, l] for l in 1:nbReturnsSegments))
	@constraint(lotSizingModel, [k in 1:nbProducts, t in 1:(tBreak - 1)],
		expBacklog[k, t] == folfOrigin[k, t] + sum(folfSlope[k, t, l] * wSegment[k, t, l] for l in 1:nbSegments) - sum(returnFolfSlope[k, l] * omegaSegment[k, l] for l in 1:nbReturnsSegments))
	## Scenario inventory
	# First, inventory evolves according to production
	@constraint(lotSizingModel, [k in 1:nbProducts, s in 1:nbScenarios],
		invRecourse[k, 1, s] - backorderRecourse[k, 1, s] == initInv[k, 1]  +  production[k, 1] - returnVol[k] - demandScenarios[k, 1, s])
	@constraint(lotSizingModel, [k in 1:nbProducts, t in 2:(tBreak-1), s in 1:nbScenarios],
		invRecourse[k, t, s] - backorderRecourse[k, t, s] == invRecourse[k, t-1, s] - backorderRecourse[k, t-1, s]  +  production[k, t] - demandScenarios[k, t, s])
	# Then inventory evolves according to recourse production
	@constraint(lotSizingModel, [k in 1:nbProducts, t in tBreak:T, s in 1:nbScenarios],
		invRecourse[k, t, s] - backorderRecourse[k, t, s] == invRecourse[k, t-1, s] - backorderRecourse[k, t-1, s] +  productionRecourse[k, t - tBreak + 1, s] - demandScenarios[k, t, s])

	# Multi-stage structure: production recourse set equal
	for t in 1:(T - tBreak + 1)
		# PRODUCTION
        nbNodes = nbProductionNodes[t]
		lengthBranch = Int(nbScenarios / nbNodes)
        for node in 1:nbNodes
			firstNodeIndex = (node-1) * lengthBranch + 1
			lastNodeIndex = node * lengthBranch
			@constraint(lotSizingModel, [k in 1:nbProducts, s in firstNodeIndex:lastNodeIndex], productionRecourse[k, t, s] == productionRecourse[k, t, firstNodeIndex])
        end
	end

	# Capacity constraint
	@constraint(lotSizingModel, [t in 1:(tBreak - 1)],
		sum(production[k, t]  for k in 1:nbProducts) <= capacity)
	@constraint(lotSizingModel, [t in 1:(T - tBreak + 1), s in 1:nbScenarios],
		sum(productionRecourse[k, t, s]  for k in 1:nbProducts) <= capacity)
	# Setup constraint
	@constraint(lotSizingModel, [k in 1:nbProducts, t in 1:(tBreak - 1)],
		production[k, t] <= capacity * setup[k, t])
	@constraint(lotSizingModel, [k in 1:nbProducts, t in 1:(T - tBreak + 1), s in 1:nbScenarios],
		productionRecourse[k, t, s] <= capacity * setup[k, t + tBreak - 1])

	# Assign wsegments to production volumes
	@constraint(lotSizingModel, [k in 1:nbProducts, t in 2:(tBreak - 1)],
		sum(wSegment[k, t, l] - wSegment[k, t-1, l] for l in 1:nbSegments) == production[k, t])
	@constraint(lotSizingModel, [k in 1:nbProducts],
		sum(wSegment[k, 1, l] for l in 1:nbSegments) == production[k, 1])
	# Logic constraints of w variables
	@constraint(lotSizingModel, [k in 1:nbProducts, t in 1:(tBreak - 1), l in 2:nbSegments],
		wSegment[k, t, l-1] >= (uRange[k, t, l] - uRange[k, t, l-1]) * lambda[k, t, l])
	@constraint(lotSizingModel, [k in 1:nbProducts, t in 1:(tBreak - 1), l in 1:nbSegments],
		wSegment[k, t, l] <= (uRange[k, t, l+1] - uRange[k, t, l]) * lambda[k, t, l])

	# Add free returns constraints
	@constraint(lotSizingModel, [k in 1:nbProducts],
		sum(omegaSegment[k, l] for l in 1:nbReturnsSegments) == returnVol[k])
	@constraint(lotSizingModel, [k in 1:nbProducts, l in 2:nbReturnsSegments],
		omegaSegment[k, l-1] >= - (vRange[k, l] - vRange[k, l-1]) * nu[k, l])
	@constraint(lotSizingModel, [k in 1:nbProducts, l in 1:nbReturnsSegments],
		omegaSegment[k, l] <= - (vRange[k, l+1] - vRange[k, l]) * nu[k, l])
	# Allow returns only if no production:
	@constraint(lotSizingModel, [k in 1:nbProducts],
		returnVol[k] <= 2 * capacity * (1-setup[k, 1]))


	# Solve optimisation model
	optimize!(lotSizingModel)

	## Return optimal decisions and expected cost
	prodOut = value.(production)
	productionRecourseOut = value.(productionRecourse)
	solveTime = MOI.get(lotSizingModel, MOI.SolveTime())

	return prodOut, productionRecourseOut, solveTime
end
