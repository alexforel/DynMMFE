function deterministic_equivalent_optim(initInv, demandForecast, capacity, nbProducts, gurobiEnv)
	## deterministic_equivalent_optim(initInv, demandForecast, capacity, nbProducts, gurobiEnv)
	# Formulates and solve a determinsitic lot-sizing problem.
	# Input:
	#   - initInv: Matrix[nbProducts], initial inventory
	#   - demandForecast: Matrix[nbProducts, T], point-estimate demand forecast
	#   - capacity: scalar, capacity limit in each period
	# Output:
	#   - prodOut: Matrix[nbProducts, T], optimal production volumes
	#   - solveTime: scalar, solve time in seconds

	if VERBOSE == 1
		println("---- Solving deterministic model ----")
	end
	## Define Optimisation model
	lotSizingModel = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(gurobiEnv), "TimeLimit" => 100, "MIPGap" => 0.01, "OutputFlag" => VERBOSE))

	## ---- Variables ----
	@variables(lotSizingModel,
		begin
			production[1:nbProducts, 1:T] >= 0
			setup[1:nbProducts, 1:T], Bin
			# Auxiliary variables
			inventory[1:nbProducts, 1:T] >= 0
			backlog[1:nbProducts, 1:T] >= 0
		end)

	## ---- Objective ----
	# Objective: min production cost, inventory costs and shortfall costs
	@objective(lotSizingModel, Min, sum(invCost' * inventory + 1000 * invCost' * backlog) + setupCost * sum(setup))

	## ---- Constraints ----
	# Calculate inventory and backlog
	@constraint(lotSizingModel, [k in 1:nbProducts, t in 2:T],
		inventory[k, t] - backlog[k, t] == production[k, t] + inventory[k, t-1] - backlog[k, t-1] - demandForecast[k, t])
	@constraint(lotSizingModel, [k in 1:nbProducts],
		inventory[k, 1] - backlog[k, 1] == production[k, 1] + initInv[k, 1]  - demandForecast[k, 1])

	# Capacity constraint
	@constraint(lotSizingModel, [t in 1:T],
		sum(production[k, t] for k in 1:nbProducts) <= capacity)
	# Setup constraint
	@constraint(lotSizingModel, [k in 1:nbProducts, t in 1:T],
		production[k, t] <= capacity * setup[k, t])

	# Solve the problem
	optimize!(lotSizingModel)
	# Return optimal values
	prodOut = value.(production)
	solveTime = MOI.get(lotSizingModel, MOI.SolveTime())

	return prodOut, solveTime
end
