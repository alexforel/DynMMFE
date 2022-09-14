## -- Sensitivity analaysis of MMFE model and value of recourse with varying capacity ---
# Reproduce results shown in EC5.1

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
using HypothesisTests

gurobiEnv = Gurobi.Env()

# Load functions
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
# Plot
include(raw"src\plot\plot_boxplot_value_of_recourse.jl")
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
NB_SIMULATIONS = 1000
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
# Scenario tree and number of scenarios
NB_SCEN_NODES = [3 6 12 24 48 48]
NB_SCENARIOS = NB_SCEN_NODES[end]
LATIN_STRATA = Matrix(undef, T, 1)
for t in 1:T
    if (t>1) && (NB_SCEN_NODES[t] == NB_SCEN_NODES[t-1])
        LATIN_STRATA[t] = LATIN_STRATA[t-1]
    else
        LATIN_STRATA[t] = latin_hypercube_strata(NB_SCEN_NODES[t])
    end
end
T_BREAK = 4
INIT_INV = 50 * ones(NB_PRODUCTS)
TIME_CORRELATION_COEFFICIENT = 0.0
PRODUCT_CORRELATION_COEFFICIENT = 0.0

## ----- Sensitivity analysis -----
VERBOSE = 0
# Sensitivity analysis parameters
capacityVector = [250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600]
uncertaintyCoeffVector = [1]
VAR_VECTOR = [400; 0.04]
forecastEvolTypeVector = ["Additive", "Multiplicative"]
demandPatternVector = ["Seasonal"]

# Create dataframe to store results
resultDataframe = DataFrame(demandPattern = String[], forecastEvolType = String[], uncertaintyCoeff =  Int64[], capacity = Int64[], model = String[], simIndex = Int64[], cost = Float64[], serviceLevel = Float64[])
start = time()
println("---------- Start of sensitivity analysis. ----------")
for demandPattern in demandPatternVector
    demandDistributions = demand_distribution_pattern(demandPattern, simulationLength)
    for forecastEvolType in forecastEvolTypeVector
        for uncertaintyCoeff in uncertaintyCoeffVector
            if forecastEvolType == "Additive"
                VAR_UPDATE = VAR_VECTOR[1, uncertaintyCoeff] * ones(T)
            elseif forecastEvolType == "Multiplicative"
                VAR_UPDATE = VAR_VECTOR[2, uncertaintyCoeff] * ones(T)
            end
            # Get true and estimated forecast evolution models
            multForecastEvolDis, multMeanVector, addMmfeMean, addMmfeCovarianceMatrix = true_forecast_evolution_and_estimated_distribution(VAR_UPDATE, demandDistributions, TIME_CORRELATION_COEFFICIENT, PRODUCT_CORRELATION_COEFFICIENT, forecastEvolType)
            if forecastEvolType == "Additive"
                trueForecastEvolDis = MvNormal(vec(addMmfeMean), addMmfeCovarianceMatrix)
            elseif forecastEvolType == "Multiplicative"
                trueForecastEvolDis = multForecastEvolDis
            end
            # Simulate forecast evolution and demand realisation
            demandRealStored, forecastStored = generate_demand_and_forecasts_from_mmfe(trueForecastEvolDis, demandDistributions, T, NB_SIMULATIONS, forecastEvolType)
            for capacity in capacityVector
                ## ------ Rolling horizon simulation ------ ##
                println("Demand pattern: $demandPattern, evolution: $forecastEvolType, uncertainty: $uncertaintyCoeff, capacity = $capacity")
                # Initialise simulation instances
                mmfeSimInstances = init_simulation_instances(NB_SIMULATIONS, INIT_INV, simulationLength, NB_PRODUCTS, T)
                extendMmfeSimInstances = init_simulation_instances(NB_SIMULATIONS, INIT_INV, simulationLength, NB_PRODUCTS, T)
                for simIndex in 1:NB_SIMULATIONS
                    print_message_iteration(simIndex)
                    for t in 1:simulationLength
                        # Solve problem using PLA
                        if forecastEvolType == "Additive"
                            mmfeProdOut, mmfePlaSolveTime = additive_mmfe_with_pla(t, forecastStored[simIndex][:, :, t], mmfeSimInstances[simIndex].inventoryState[:, t], capacity, NB_PRODUCTS, addMmfeMean, addMmfeCovarianceMatrix)
                        elseif forecastEvolType == "Multiplicative"
                            mmfeProdOut, mmfePlaSolveTime = multiplicative_mmfe_with_pla(t, forecastStored[simIndex][:, :, t], mmfeSimInstances[simIndex].inventoryState[:, t], multMeanVector, scale(multForecastEvolDis), capacity, NB_PRODUCTS)
                        end
                        extendPlaProdOut, extendProdRecourse, extendPlaSolveTime = extended_mmfe_with_pla(t, forecastStored[simIndex][:, :, t], extendMmfeSimInstances[simIndex].inventoryState[:, t], forecastEvolType, trueForecastEvolDis, capacity, T_BREAK, NB_PRODUCTS)
                        # Update simulation instances
                        update_instance_state!(t, mmfeProdOut, mmfePlaSolveTime, demandRealStored[simIndex][:, t], mmfeSimInstances[simIndex], NB_PRODUCTS)
                        update_instance_state!(t, extendPlaProdOut, extendProdRecourse, extendPlaSolveTime, demandRealStored[simIndex][:, t], extendMmfeSimInstances[simIndex], NB_PRODUCTS)
                    end
                end
                # KPI calculation
                # Calculate realised cost over simulation
                mmfePlaCosts, addMmfeSl, addMmfeAggSl = kpis(mmfeSimInstances, backlogCost, demandRealStored, NB_PRODUCTS, NB_SIMULATIONS)
                extendMmfePlaCosts, extendMmfeSl, extendMmfeAggSl = kpis(extendMmfeSimInstances, backlogCost, demandRealStored, NB_PRODUCTS, NB_SIMULATIONS)
                for s in 1:NB_SIMULATIONS
                    push!(resultDataframe, [demandPattern forecastEvolType uncertaintyCoeff capacity "mmfe" s mmfePlaCosts[s, 4] addMmfeAggSl[s]])
                    push!(resultDataframe, [demandPattern forecastEvolType uncertaintyCoeff capacity "extended" s extendMmfePlaCosts[s, 4] extendMmfeAggSl[s]])
                end
            end
        end
    end
end
elapsed = time() - start
println("\n---------- End of sensitivity analysis simulations, elapsed: $elapsed s. ----------")
# Export dataframe result to csv
CSV.write("output\\resultDataframe_capacity_sensitivity.csv", resultDataframe)

## ---- Calculate average cost and value of recourse ----
addResultDf = filter(row -> row[:forecastEvolType] == "Additive", resultDataframe)
multResultDf = filter(row -> row[:forecastEvolType] == "Multiplicative", resultDataframe)
# Cost of mmfe and extended model as a function of available capacity
mmfeAvgCost = zeros(length(capacityVector))
extMmfeAvgCost = zeros(length(capacityVector))
mmfeConf = zeros(length(capacityVector)) # store confidence interval on average cost
extMmfeConf = zeros(length(capacityVector)) # store confidence interval on average cost
relativeCostSamples = zeros(NB_SIMULATIONS, length(capacityVector))
relativeCostConf = zeros(length(capacityVector)) # store confidence interval on relative cost
for c in 1:length(capacityVector)
    cap = capacityVector[c]
    # Filter dataframe to relevant instance
    tempDf = filter(row -> row[:capacity] == cap, multResultDf)
    tempMmfeDf = filter(row -> row[:model] == "mmfe", tempDf)
    # Get costs
    mmfeSamples = tempMmfeDf[!, :cost]
    mmfeAvgCost[c] = mean(mmfeSamples)
    mmfeConf[c] = 1.96 * std(mmfeSamples) / sqrt(length(mmfeSamples))
    tempExtMmfeDf = filter(row -> row[:model] == "extended", tempDf)
    extSamples = tempExtMmfeDf[!, :cost]
    extMmfeAvgCost[c] = mean(extSamples)
    extMmfeConf[c] = 1.96 * std(extSamples) / sqrt(length(extSamples))
    # Calculate relative cost
    relativeCostSamples[:, c] = extSamples ./ mmfeSamples * 100
    println("c = ", c, "pvalue = ", pvalue(OneSampleTTest(relativeCostSamples[:, c], 100)))
    relativeCostConf[c] = 1.96 * std(relativeCostSamples[:, c]) / sqrt(length(relativeCostSamples[:, c]))
end
# Visualise results
plot(capacityVector, [mmfeAvgCost extMmfeAvgCost], yerror = [mmfeConf extMmfeConf])
relativeCost = mean(relativeCostSamples, dims = 1)
plot(capacityVector, relativeCost', yerror = relativeCostConf)
# Store results in csv file
capaSensDf = DataFrame([capacityVector mmfeAvgCost mmfeConf extMmfeAvgCost extMmfeConf relativeCost' relativeCostConf])
CSV.write("output\\capaSensAdditive.csv", capaSensDf)
