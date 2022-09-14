function pla_lot_sizing_free_returns_no_prod_recourse(initInv, uRange, folfOrigin, folfSlope, complFolfOrigin, complFolfSlope, vRange, returnFolfSlope, returnComplFolfSlope, nbSegments, nbReturnsSegments, backlogCost, capacity, nbProducts, gurobiEnv)
	## pla_lot_sizing_no_prod_recourse(initInv, uRange, folfOrigin, folfSlope, complFolfOrigin, complFolfSlope, nbSegments, backlogCost, capacity, nbProducts, gurobiEnv)
	# Formulate a stochastic lot-sizing with piecewise-linearisation of the inventory and backlog functions and solve it using Gurobi.
	# Input:
	#   - initInv, uRange, folfOrigin, folfSlope, complFolfOrigin, complFolfSlope, nbSegments, backlogCost, capacity, nbProducts, gurobiEnv
	# Output:
	#   - prodOut, solveTime

	if VERBOSE == 1
		println("---- Solving PLA model ----")
	end
	## Define Optimisation model
	lotSizingModel = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(gurobiEnv), "TimeLimit" => 100, "MIPGap" => 0.01, "OutputFlag" => VERBOSE))

	## ---- Variables ----
	@variables(lotSizingModel,
		begin
			production[1:nbProducts, 1:T] >= 0
			setup[1:nbProducts, 1:T], Bin
			# Auxiliary variables
			lambda[1:nbProducts, 1:T, 1:nbSegments], Bin
			expInv[1:nbProducts, 1:T]
			expBacklog[1:nbProducts, 1:T]
			wSegment[1:nbProducts, 1:T, 1:nbSegments] >= 0
			## Free return variables
			returnVol[1:nbProducts] >= 0
			# Auxiliary variables
			nu[1:nbProducts, 1:nbReturnsSegments], Bin
			omegaSegment[1:nbProducts, 1:nbReturnsSegments] >= 0
		end)

	## ---- Objective ----
	# Objective: min production cost, inventory costs and shortfall costs
	@objective(lotSizingModel, Min, sum(sum(invCost[k] * expInv[k,t] + backlogCost[k] * expBacklog[k,t] + setupCost * setup[k,t] for t in 1:T) for k in 1:nbProducts))

	## ---- Constraints ----
	## PLA-specific constraints
	# Expected inventory and backlog
	@constraint(lotSizingModel, [k in 1:nbProducts, t in 1:T],
		expInv[k, t] == complFolfOrigin[k, t] + sum(complFolfSlope[k, t, l] * wSegment[k, t, l] for l in 1:nbSegments) - sum(returnComplFolfSlope[k, l] * omegaSegment[k, l] for l in 1:nbReturnsSegments))
	@constraint(lotSizingModel, [k in 1:nbProducts, t in 1:T],
		expBacklog[k, t] == folfOrigin[k, t] + sum(folfSlope[k, t, l] * wSegment[k, t, l] for l in 1:nbSegments) - sum(returnFolfSlope[k, l] * omegaSegment[k, l] for l in 1:nbReturnsSegments))
	# Assign wSegments to production volumes
	@constraint(lotSizingModel, [k in 1:nbProducts, t in 2:T],
		sum(wSegment[k, t, l] - wSegment[k, t-1, l] for l in 1:nbSegments) == production[k, t])
	@constraint(lotSizingModel, [k in 1:nbProducts],
		sum(wSegment[k, 1, l] for l in 1:nbSegments) == production[k, 1])
	# Logic constraints of w variables
	@constraint(lotSizingModel, [k in 1:nbProducts, t in 1:T, l in 2:nbSegments],
		wSegment[k, t, l-1] >= (uRange[k, t, l] - uRange[k, t, l-1]) * lambda[k, t, l])
	@constraint(lotSizingModel, [k in 1:nbProducts, t in 1:T, l in 1:nbSegments],
		wSegment[k, t, l] <= (uRange[k, t, l+1] - uRange[k, t, l]) * lambda[k, t, l])
	## Lot-sizing specific constraints
	# Capacity constraint
	@constraint(lotSizingModel, [t in 1:T],
		sum(production[k, t] for k in 1:nbProducts) <= capacity)
	# Setup constraint
	@constraint(lotSizingModel, [k in 1:nbProducts, t in 1:T],
		production[k, t] <= capacity * setup[k, t])

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

	# Solve problem
	optimize!(lotSizingModel)

	# Read and output solution of variables
	prodOut = value.(production)
	solveTime = MOI.get(lotSizingModel, MOI.SolveTime())

	return prodOut, solveTime
end
