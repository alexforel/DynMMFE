## -- Sensitivity analaysis of scenario structure. --
#  Reproduce results provided in EC5.2

## Initialisation
# Set working directory
cd("C:/Users/Alexandre//Documents/Code/DynMMFE")

using LinearAlgebra
using Random, Distributions, Statistics
using Plots.PlotMeasures
using StatsPlots
using QuadGK
using JuMP, Gurobi
using DataFrames
using CSV
using LaTeXTabulars

gurobiEnv = Gurobi.Env()

## Load functions
# Types
include(raw"src\types\PiecewiseLinear.jl")
include(raw"src\types\SimulationInstance.jl")
include(raw"src\types\init_simulation_instances.jl")
# First-order loss function
include(raw"src\folf\folf_linearisation.jl")
include(raw"src\folf\first_order_loss_function.jl")
# Add-MMFE
include(raw"src\mmfe\additive_mmfe_with_pla.jl")
include(raw"src\mmfe\cumulative_demand_distribution_matrix_from_additive_mmfe.jl")
# Mult-MMFE
include(raw"src\mmfe\multiplicative_mmfe_with_pla.jl")
include(raw"src\mmfe\cumulative_demand_distribution_matrix_from_multiplicative_mmfe.jl")
include(raw"src\mmfe\fenton_wilkinson_approx_sum_log_normal.jl")
# Extended MMFE
include(raw"src\mmfe\extended_mmfe_with_pla.jl")
include(raw"src\mmfe\build_scenario_tree_from_forecast_evolution.jl")
# Methods
include(raw"src\methods\run_single_period_optimisation_and_update!.jl")
include(raw"src\methods\update_demand_and_forecast!.jl")
include(raw"src\methods\update_instance_state!.jl")
# Optimisation problems
include(raw"src\optim\pla_lot_sizing_no_prod_recourse.jl")
include(raw"src\optim\deterministic_equivalent_optim.jl")
include(raw"src\optim\extended_pla_scenario_lot_sizing_with_prod_recourse.jl")
# KPIs
include(raw"src\kpi\kpis.jl")
include(raw"src\kpi\realised_costs.jl")
include(raw"src\kpi\realised_service_level.jl")
# Distribution
include(raw"src\distribution\true_forecast_evolution_and_estimated_distribution.jl")
# Demand
include(raw"src\demand\demand_distribution_pattern.jl")
include(raw"src\demand\generate_demand_and_forecasts_from_mmfe.jl")
# Output
include(raw"src\output\filter_sensitivity_results.jl")
# Misc
include(raw"src\print_message_iteration.jl")
include(raw"src\simpson_integration.jl")
include(raw"src\sampling\latin_hypercube_sampling_with_multidimensional_uniformity.jl")

# Piecewise-linear approximation of the inverse of the standard folf
xLim = 5
first_order_loss_function_standard_normal(x) = quadgk(t -> (1 - cdf(Normal(), t)), x, Inf)[1]
epsilonRange = collect(range(-xLim, stop = xLim, length = 2000))
pwlaStandardFolf = PiecewiseLinear(epsilonRange, first_order_loss_function_standard_normal)

## Simulation environment parameters
Random.seed!(41193) # Fix random seed
NB_PRODUCTS = 2
T = 6
# Simulation length
preSimT = 0 # warm-up before demand
inSimT = 12 # periods with demand
postSimT = 0 # cool-down after demand
simulationLength = preSimT + inSimT + postSimT
# Cost parameters
invCost = 1 .+ 0.5 * rand(NB_PRODUCTS) # ones(NB_PRODUCTS)
backlogCost = 10 * invCost
setupCost = 150
# Model parameters
NB_SEGMENTS = 40 # should be an even number
INIT_INV = 50 * ones(NB_PRODUCTS)
TIME_CORRELATION_COEFFICIENT = 0.0
PRODUCT_CORRELATION_COEFFICIENT = 0.0

## ----- Sensitivity analysis -----
NB_SIMULATIONS = 1000
VERBOSE = 0
# Sensitivity analysis parameters
capacity = 300
VAR_VECTOR = [400, 0.04]
demandPattern = "Seasonal"
demandDistributions = demand_distribution_pattern(demandPattern, simulationLength)
# Scenario tree and number of scenarios
scenarioTreeVector = [[2 4 8 8 16 16], [3 6 12 24 48 48], [3 9 18 36 36 72]]
forecastEvolVector = ["Additive", "Multiplicative"]
tBreakVector = [3, 4, 5, 6, 7]

# Create dataframe to store results
resultDataframe = DataFrame(treeIndex = Int64[], forecastEvolType = String[], tbreak = Int64[], simIndex = Int64[], cost = Float64[], serviceLevel = Float64[], solveTime = Float64[])
start = time()
for treeIndex in 1:length(scenarioTreeVector)
    println("treeIndex = $treeIndex")
    global NB_SCEN_NODES = scenarioTreeVector[treeIndex]
    global NB_SCENARIOS = NB_SCEN_NODES[end]
    global LATIN_STRATA = Matrix(undef, T, 1)
    for t in 1:T
        if (t>1) && (NB_SCEN_NODES[t] == NB_SCEN_NODES[t-1])
            LATIN_STRATA[t] = LATIN_STRATA[t-1]
        else
            LATIN_STRATA[t] = latin_hypercube_strata(NB_SCEN_NODES[t])
        end
    end
    for forecastEvolType in forecastEvolVector
        println("forecastEvolType = $forecastEvolType")
        # Get true forecast evolution distribution
        if forecastEvolType == "Additive"
            VAR_UPDATE = VAR_VECTOR[1, 1] * ones(T)
        elseif forecastEvolType == "Multiplicative"
            VAR_UPDATE = VAR_VECTOR[2, 1] * ones(T)
        end
        multForecastEvolDis, multMeanVector, addMmfeMean, addMmfeCovarianceMatrix = true_forecast_evolution_and_estimated_distribution(VAR_UPDATE, demandDistributions, TIME_CORRELATION_COEFFICIENT, PRODUCT_CORRELATION_COEFFICIENT, forecastEvolType, 1)
        if forecastEvolType == "Additive"
            trueForecastEvolDis = MvNormal(vec(addMmfeMean), addMmfeCovarianceMatrix)
        elseif forecastEvolType == "Multiplicative"
            trueForecastEvolDis = multForecastEvolDis
        end
        # Simulate forecast evolution and demand realisation
        demandRealStored, forecastStored = generate_demand_and_forecasts_from_mmfe(trueForecastEvolDis, demandDistributions, T, NB_SIMULATIONS, forecastEvolType)
        for tBreak in tBreakVector
            println("tBreak+1 = $tBreak")
            # Initialise simulation instances
            extendMmfeSimInstances = init_simulation_instances(NB_SIMULATIONS, INIT_INV, simulationLength, NB_PRODUCTS, T)
            ## ------ Rolling horizon simulation ------ ##
            for simIndex in 1:NB_SIMULATIONS
                print_message_iteration(simIndex)
                for t in 1:simulationLength
                    if tBreak == 7
                        if forecastEvolType == "Additive"
                            extendPlaProdOut, extendPlaSolveTime = additive_mmfe_with_pla(t, forecastStored[simIndex][:, :, t], extendMmfeSimInstances[simIndex].inventoryState[:, t], capacity, NB_PRODUCTS, addMmfeMean, addMmfeCovarianceMatrix)
                        elseif forecastEvolType == "Multiplicative"
                            extendPlaProdOut, extendPlaSolveTime = multiplicative_mmfe_with_pla(t, forecastStored[simIndex][:, :, t], extendMmfeSimInstances[simIndex].inventoryState[:, t], multMeanVector, scale(multForecastEvolDis), capacity, NB_PRODUCTS)
                        end
                        update_instance_state!(t, extendPlaProdOut, extendPlaSolveTime, demandRealStored[simIndex][:, t], extendMmfeSimInstances[simIndex], NB_PRODUCTS)
                    else
                        extendPlaProdOut, extendProdRecourse, extendPlaSolveTime = extended_mmfe_with_pla(t, forecastStored[simIndex][:, :, t], extendMmfeSimInstances[simIndex].inventoryState[:, t], forecastEvolType, trueForecastEvolDis, capacity, tBreak, NB_PRODUCTS)
                        update_instance_state!(t, extendPlaProdOut, extendProdRecourse, extendPlaSolveTime, demandRealStored[simIndex][:, t], extendMmfeSimInstances[simIndex], NB_PRODUCTS)
                    end
                end
                # KPI calculation
                extendMmfePlaCosts, extendMmfeSl, extendMmfeAggSl = kpis(extendMmfeSimInstances[simIndex], backlogCost, demandRealStored[simIndex], NB_PRODUCTS)
                tempSolveTime = sum(extendMmfeSimInstances[simIndex].solveTime)
                # Store result in Dataframe
                push!(resultDataframe, [treeIndex forecastEvolType tBreak simIndex extendMmfePlaCosts[4] extendMmfeAggSl tempSolveTime])
            end
        end
    end
end
elapsed = time() - start
println("\n---------- End of sensitivity analysis simulations, elapsed: $elapsed s. ----------")
# Export dataframe result to csv
CSV.write("output\\sensitivityScenarioTreeTbreak.csv", resultDataframe)

## Analyse results
addresultDf = filter(row -> row[:forecastEvolType] == "Additive", resultDataframe)
multresultDf = filter(row -> row[:forecastEvolType] == "Multiplicative", resultDataframe)
# Calculate performance indicators: value of recourse and solution time
addaverageRelativeCost, addcostConfInterval, addsolveTimeVector, addtimeConfInterval = filter_sensitivity_results(addresultDf)
multaverageRelativeCost, multcostConfInterval, multsolveTimeVector, multtimeConfInterval = filter_sensitivity_results(multresultDf)
# Plot results
plot(tBreakVector, addaverageRelativeCost', yerror = addcostConfInterval')
plot(tBreakVector, addsolveTimeVector', yerror = addtimeConfInterval')
plot(tBreakVector, multaverageRelativeCost', yerror = multcostConfInterval')
plot(tBreakVector, multsolveTimeVector', yerror = multtimeConfInterval')
# Export results to csv table
tBreakSensitivityTable = DataFrame([tBreakVector.-1 addaverageRelativeCost' addcostConfInterval' addsolveTimeVector' addtimeConfInterval'], :auto)
CSV.write("output\\csv\\additiveTreeSensitivity.csv", tBreakSensitivityTable)
tBreakSensitivityTable = DataFrame([tBreakVector.-1 multaverageRelativeCost' multcostConfInterval' multsolveTimeVector' multtimeConfInterval'], :auto)
CSV.write("output\\csv\\multiplicativeTreeSensitivity.csv", tBreakSensitivityTable)
