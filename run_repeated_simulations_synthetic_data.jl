## ---- Main -----
# Run rolling-horizon simulations for artifical demand and forecast data. Forecast evolution follows the Martingale model of forecast evolution. We compare the performances of the right MMFE model and a MMFE model with wrong distribution assumptions. The true forecast evolution process is set to follow either the additive or multiplicative MMFE.
# Reproduce results presented in Section 5.1.2, Section 5.1.3, EC.3 and EC.4

# Structure:
# - Define true forecast evolution process
# - Estimate parameters of mismatched MMFE model from simulations of true MMFE process
# - Create demand patterns over the simulation length
# - Perfom rolling-horizon simulations
# - Visualise and export results and performance indicators

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
# Sampling
include(raw"src\sampling\latin_hypercube_sampling_with_multidimensional_uniformity.jl")
# Plot
include(raw"src\plot\plot_mean_cost_bar_over_simulation_enviroments.jl")
include(raw"src\plot\plot_boxplot_value_of_recourse.jl")
include(raw"src\plot\plot_boxplot_value_of_recourse_over_simulation_enviroments.jl")
include(raw"src\plot\plot_mean_bar_value_of_recourse_over_simulation_enviroments.jl")
# Output
include(raw"src\output\average_and_relative_costs_per_simulation_environment.jl")
include(raw"src\output\export_average_results_to_latex_table.jl")
include(raw"src\output\export_value_of_recourse_to_latex_table.jl")
include(raw"src\output\interaction_effect_capacity_uncertainty.jl")
# Misc
include(raw"src\print_message_iteration.jl")
include(raw"src\simpson_integration.jl")
include(raw"src\demand_driven_stochastic.jl")

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
simulationLength = 12 # periods with demand
# Cost parameters
invCost = 1 .+ 0.5 * rand(NB_PRODUCTS) # ones(NB_PRODUCTS)
backlogCost = 10 * invCost
setupCost = 150
# Model parameters
NB_SEGMENTS = 40 # should be an even number
# Scenario tree and number of scenarios
NB_SCEN_NODES = [3 6 12 24 48 48] # number of nodes in each period of the horizon
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

## ----- Rolling-horizon simulations -----
VERBOSE = 0
# Sensitivity analysis parameters
capacityVector = [300 ; 500]
uncertaintyCoeffVector = [1, 2, 3]
VAR_VECTOR = [100 400 700 ; 0.01 0.04 0.07]
forecastEvolTypeVector = ["Additive", "Multiplicative"]
demandPatternVector = ["Stationary", "Random", "Seasonal"]

# Create dataframe to store results
resultDataframe = DataFrame(demandPattern = String[], forecastEvolType = String[], uncertaintyCoeff =  Int64[], capacity = Int64[], model = String[], simIndex = Int64[], cost = Float64[], serviceLevel = Float64[])
start = time()
println("---------- Start of sensitivity analysis. ----------")
for demandPattern in demandPatternVector
    # Demand pattern and forecast evolution model: true model and estimated model
    demandDistributions = demand_distribution_pattern(demandPattern, simulationLength)
    for forecastEvolType in forecastEvolTypeVector
        for uncertaintyCoeff in uncertaintyCoeffVector
            if forecastEvolType == "Additive"
                VAR_UPDATE = VAR_VECTOR[1, uncertaintyCoeff] * ones(T)
            elseif forecastEvolType == "Multiplicative"
                VAR_UPDATE = VAR_VECTOR[2, uncertaintyCoeff] * ones(T)
            end
            multForecastEvolDis, multMeanVector, addMmfeMean, addMmfeCovarianceMatrix, ddsMean, ddsStd = true_forecast_evolution_and_estimated_distribution(VAR_UPDATE, demandDistributions, TIME_CORRELATION_COEFFICIENT, PRODUCT_CORRELATION_COEFFICIENT, forecastEvolType)
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
                detSimInstances = init_simulation_instances(NB_SIMULATIONS, INIT_INV, simulationLength, NB_PRODUCTS, T)
                ddStochSimInstances = init_simulation_instances(NB_SIMULATIONS, INIT_INV, simulationLength, NB_PRODUCTS, T)
                multMmfeSimInstances = init_simulation_instances(NB_SIMULATIONS, INIT_INV, simulationLength, NB_PRODUCTS, T)
                addMmfeSimInstances = init_simulation_instances(NB_SIMULATIONS, INIT_INV, simulationLength, NB_PRODUCTS, T)
                extendMmfeSimInstances = init_simulation_instances(NB_SIMULATIONS, INIT_INV, simulationLength, NB_PRODUCTS, T)
                for simIndex in 1:NB_SIMULATIONS
                    print_message_iteration(simIndex)
                    for t in 1:simulationLength
                        run_single_period_optimisation_and_update!(t, forecastStored[simIndex][:, :, t], demandRealStored[simIndex][:, t], detSimInstances[simIndex], addMmfeSimInstances[simIndex], addMmfeMean, addMmfeCovarianceMatrix, multMmfeSimInstances[simIndex], multForecastEvolDis, multMeanVector, extendMmfeSimInstances[simIndex], forecastEvolType, trueForecastEvolDis, capacity, T_BREAK, NB_PRODUCTS)
                        # Stochastic benchmark - Demand uncertainty
                        ddsForecast = zeros(NB_PRODUCTS, T)
                        ddsStdVec = zeros(NB_PRODUCTS, T)
                        for tau in 1:T
                            indexToRead = t+tau-1
                            if indexToRead <= simulationLength
                                ddsForecast[:, tau] = ddsMean[:, indexToRead]
                                ddsStdVec[:, tau] = ddsStd[:, indexToRead]
                            end
                        end
                        ddsProdOut, ddsSolveTime = demand_driven_stochastic(t, ddsForecast, ddsStdVec, ddStochSimInstances[simIndex].inventoryState[:, t], capacity, NB_PRODUCTS, forecastEvolType)
                        update_instance_state!(t, ddsProdOut, ddsSolveTime, demandRealStored[simIndex][:, t], ddStochSimInstances[simIndex], NB_PRODUCTS)
                    end
                end
                ## KPI calculation
                # Calculate realised cost over simulation
                detCosts, _, detAggSl = kpis(detSimInstances, backlogCost, demandRealStored, NB_PRODUCTS, NB_SIMULATIONS)
                ddsCosts, _, ddsAggSl = kpis(ddStochSimInstances, backlogCost, demandRealStored, NB_PRODUCTS, NB_SIMULATIONS)
                addMmfePlaCosts, _, addMmfeAggSl = kpis(addMmfeSimInstances, backlogCost, demandRealStored, NB_PRODUCTS, NB_SIMULATIONS)
                multMmfePlaCosts, _, multMmfeAggSl = kpis(multMmfeSimInstances, backlogCost, demandRealStored, NB_PRODUCTS, NB_SIMULATIONS)
                extendMmfePlaCosts, _, extendMmfeAggSl = kpis(extendMmfeSimInstances, backlogCost, demandRealStored, NB_PRODUCTS, NB_SIMULATIONS)
                # Store result in Dataframe
                for s in 1:NB_SIMULATIONS
                    push!(resultDataframe, [demandPattern forecastEvolType uncertaintyCoeff capacity "deterministic" s detCosts[s, 4] detAggSl[s]])
                    push!(resultDataframe, [demandPattern forecastEvolType uncertaintyCoeff capacity "dd-stoch" s ddsCosts[s, 4] ddsAggSl[s]])
                    push!(resultDataframe, [demandPattern forecastEvolType uncertaintyCoeff capacity "additive" s addMmfePlaCosts[s, 4] addMmfeAggSl[s]])
                    push!(resultDataframe, [demandPattern forecastEvolType uncertaintyCoeff capacity "multiplicative" s multMmfePlaCosts[s, 4] multMmfeAggSl[s]])
                    push!(resultDataframe, [demandPattern forecastEvolType uncertaintyCoeff capacity "extended" s extendMmfePlaCosts[s, 4] extendMmfeAggSl[s]])
                end
            end
        end
    end
end
elapsed = time() - start
println("\n---------- End of sensitivity analysis simulations, elapsed: $elapsed s. ----------")
# Export dataframe result to csv
CSV.write("output\\resultDataframe_synth_with_benchmark.csv", resultDataframe)
# resultDataframe = CSV.read("output\\resultDataframe_synth_with_benchmark.csv", DataFrame)

## ---- Results analysis ----
addResultDf = filter(row -> row[:forecastEvolType] == "Additive", resultDataframe)
multResultDf = filter(row -> row[:forecastEvolType] == "Multiplicative", resultDataframe)
# Overall boxplot
@df addResultDf boxplot(:model, :cost)
@df multResultDf boxplot(:model, :cost)
@df addResultDf boxplot(:model, :serviceLevel)
@df multResultDf boxplot(:model, :serviceLevel)

# Boxplot of additive and multiplicative over all simulation environments
plot_mean_cost_bar_over_simulation_enviroments(addResultDf)
plot_mean_cost_bar_over_simulation_enviroments(multResultDf)

# Export average results in CSV with a table format
NB_MODELS = 5
addAvgCost, addRelCost, addStatSignRelativeCost = average_and_relative_costs_per_simulation_environment(addResultDf)
multAvgCost, multRelCost, multStatSignRelativeCost = average_and_relative_costs_per_simulation_environment(multResultDf)
export_average_results_to_latex_table(addAvgCost, addRelCost, addStatSignRelativeCost, multAvgCost, multRelCost, multStatSignRelativeCost)

## ---- Value of recourse ----
# plot_boxplot_value_of_recourse(addResultDf, multResultDf)
plot_boxplot_value_of_recourse_over_simulation_enviroments(addResultDf, "additive")
plot_boxplot_value_of_recourse_over_simulation_enviroments(multResultDf, "multiplicative")
plot_mean_bar_value_of_recourse_over_simulation_enviroments(addResultDf, "additive")
plot_mean_bar_value_of_recourse_over_simulation_enviroments(multResultDf, "multiplicative")
# Export average and quartile results to CSV
export_value_of_recourse_to_latex_table(addResultDf, "additive", "addValueRecourseTable")
export_value_of_recourse_to_latex_table(multResultDf, "multiplicative", "multValueRecourseTable")

## Interaction of simulation parameters
# For each demand pattern, show the value of forecast evolution as a function of the uncertainty parameter for the two capacity settings
addAverageRelativeCost, addCostConfInterval = interaction_effect_capacity_uncertainty("additive", addResultDf)
multAverageRelativeCost, multCostConfInterval = interaction_effect_capacity_uncertainty("multiplicative", multResultDf)
plot(1:3, addAverageRelativeCost, title = "Aditive MMFE", ribbon = (addCostConfInterval, addCostConfInterval))
plot(1:3, multAverageRelativeCost, ribbon = (multCostConfInterval, multCostConfInterval), title = "Multiplicative MMFE")
# Export to csv file
addInteractionDf = DataFrame([1:3 100 .- addAverageRelativeCost addCostConfInterval])
CSV.write("output\\addInteraction.csv", addInteractionDf)
multInteractionDf = DataFrame([1:3 100 .- multAverageRelativeCost multCostConfInterval])
CSV.write("output\\multInteraction.csv", multInteractionDf)

# Extended models with recourse
addExtAverageRelativeCost, addExtCostConfInterval = interaction_effect_capacity_uncertainty("extended", addResultDf)
multExtAverageRelativeCost, multExtCostConfInterval = interaction_effect_capacity_uncertainty("extended", multResultDf)
# Export to csv file
addExtInteractionDf = DataFrame([1:3 100 .- addExtAverageRelativeCost addExtCostConfInterval])
CSV.write("output\\addExtInteraction.csv", addExtInteractionDf)
mulExttInteractionDf = DataFrame([1:3 100 .- multExtAverageRelativeCost multExtCostConfInterval])
CSV.write("output\\multExtInteraction.csv", mulExttInteractionDf)
