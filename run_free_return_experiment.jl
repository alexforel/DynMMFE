## ---- Experiment with free returns -----
# Reproduce results in EC.6

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
include(raw"src\optim\pla_lot_sizing_free_returns_no_prod_recourse.jl")
include(raw"src\optim\extended_pla_scenario_lot_sizing_with_prod_recourse.jl")
include(raw"src\optim\extended_pla_scenario_lot_sizing_free_returns_with_prod_recourse.jl")
# KPIs
include(raw"src\kpi\kpis.jl")
include(raw"src\kpi\realised_costs.jl")
include(raw"src\kpi\realised_service_level.jl")
# Distribution
include(raw"src\distribution\true_forecast_evolution_and_estimated_distribution.jl")
# Demand
include(raw"src\demand\demand_distribution_pattern.jl")
include(raw"src\demand\generate_demand_and_forecasts_from_mmfe.jl")
# Misc
include(raw"src\print_message_iteration.jl")
include(raw"src\simpson_integration.jl")
include(raw"src\sampling\latin_hypercube_sampling_with_multidimensional_uniformity.jl")
# Output
include(raw"src\output\free_returns_format_results.jl")

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

## ----- Sensitivity analysis -----
VERBOSE = 0
demandPattern = "Seasonal"
capacity = 300
uncertaintyCoeff = 1
VAR_VECTOR = [400; 0.04]
# Sensitivity analysis parameters
forecastEvolTypeVector = ["Additive", "Multiplicative"]
timeCorr = 0.0
prodCorr = 0.0

# Global parameter unique to free returns analysis
NB_RETURNS_SEGMENTS = 40

# Create dataframe to store results
resultDataframe = DataFrame(demandPattern = String[], forecastEvolType = String[],
                            model = String[], simIndex = Int64[], cost = Float64[],
                            serviceLevel = Float64[])

start = time()
# Demand pattern and forecast evolution model: true model and estimated model
demandDistributions = demand_distribution_pattern(demandPattern, simulationLength)
println("---------- Start of simulations with free return. ----------")
for forecastEvolType in forecastEvolTypeVector# Get true and estimated forecast evolution models
    println("---- MMFE evolution: $forecastEvolType ----")
    if forecastEvolType == "Additive"
        VAR_UPDATE = VAR_VECTOR[1, uncertaintyCoeff] * ones(T)
    elseif forecastEvolType == "Multiplicative"
        VAR_UPDATE = VAR_VECTOR[2, uncertaintyCoeff] * ones(T)
    end
    multForecastEvolDis, multMeanVector, addMmfeMean, addMmfeCovarianceMatrix, ddsMean, ddsStd = true_forecast_evolution_and_estimated_distribution(VAR_UPDATE, demandDistributions, timeCorr, prodCorr, forecastEvolType)
    if forecastEvolType == "Additive"
        trueForecastEvolDis = MvNormal(vec(addMmfeMean), addMmfeCovarianceMatrix)
    elseif forecastEvolType == "Multiplicative"
        trueForecastEvolDis = multForecastEvolDis
    end
    # Simulate forecast evolution and demand realisation
    demandRealStored, forecastStored = generate_demand_and_forecasts_from_mmfe(trueForecastEvolDis, demandDistributions, T, NB_SIMULATIONS, forecastEvolType)
    ## ------ Rolling horizon simulation ------ ##
    # Initialise simulation instances
    mmfeSimInstances = init_simulation_instances(NB_SIMULATIONS, INIT_INV, simulationLength, NB_PRODUCTS, T)
    extendMmfeSimInstances = init_simulation_instances(NB_SIMULATIONS, INIT_INV, simulationLength, NB_PRODUCTS, T)
    freeReturnSimInstances = init_simulation_instances(NB_SIMULATIONS, INIT_INV, simulationLength, NB_PRODUCTS, T)
    extendFreeReturnSimInstances = init_simulation_instances(NB_SIMULATIONS, INIT_INV, simulationLength, NB_PRODUCTS, T)
    misspecSimInstances = init_simulation_instances(NB_SIMULATIONS, INIT_INV, simulationLength, NB_PRODUCTS, T)
    misspecFreeReturnSimInstances = init_simulation_instances(NB_SIMULATIONS, INIT_INV, simulationLength, NB_PRODUCTS, T)
    for simIndex in 1:NB_SIMULATIONS
        print_message_iteration(simIndex)
        for t in 1:simulationLength
            # Solve problem using PLA
            if forecastEvolType == "Additive"
                mmfeProdOut, mmfePlaSolveTime = additive_mmfe_with_pla(t, forecastStored[simIndex][:, :, t], mmfeSimInstances[simIndex].inventoryState[:, t], capacity, NB_PRODUCTS, addMmfeMean, addMmfeCovarianceMatrix)
                freeReturnMmfeProdOut, freeReturnMmfePlaSolveTime = additive_mmfe_with_pla(t, forecastStored[simIndex][:, :, t], freeReturnSimInstances[simIndex].inventoryState[:, t], capacity, NB_PRODUCTS, addMmfeMean, addMmfeCovarianceMatrix, useFreeReturn=true)
                # Mis-specified MMFE model
                misspecMmfeProdOut, misspecMmfePlaSolveTime = multiplicative_mmfe_with_pla(t, forecastStored[simIndex][:, :, t], misspecSimInstances[simIndex].inventoryState[:, t], multMeanVector, scale(multForecastEvolDis), capacity, NB_PRODUCTS)
                misspecFreeReturnMmfeProdOut, misspecFreeReturnMmfePlaSolveTime = multiplicative_mmfe_with_pla(t, forecastStored[simIndex][:, :, t], misspecFreeReturnSimInstances[simIndex].inventoryState[:, t], multMeanVector, scale(multForecastEvolDis), capacity, NB_PRODUCTS, useFreeReturn=true)
            elseif forecastEvolType == "Multiplicative"
                mmfeProdOut, mmfePlaSolveTime = multiplicative_mmfe_with_pla(t, forecastStored[simIndex][:, :, t], mmfeSimInstances[simIndex].inventoryState[:, t], multMeanVector, scale(multForecastEvolDis), capacity, NB_PRODUCTS)
                freeReturnMmfeProdOut, freeReturnMmfePlaSolveTime = multiplicative_mmfe_with_pla(t, forecastStored[simIndex][:, :, t], freeReturnSimInstances[simIndex].inventoryState[:, t], multMeanVector, scale(multForecastEvolDis), capacity, NB_PRODUCTS, useFreeReturn=true)
                # Mis-specified MMFE model - Free returns but no recourse
                misspecMmfeProdOut, misspecMmfePlaSolveTime = additive_mmfe_with_pla(t, forecastStored[simIndex][:, :, t], misspecSimInstances[simIndex].inventoryState[:, t], capacity, NB_PRODUCTS, addMmfeMean, addMmfeCovarianceMatrix)
                misspecFreeReturnMmfeProdOut, misspecFreeReturnMmfePlaSolveTime = additive_mmfe_with_pla(t, forecastStored[simIndex][:, :, t], misspecFreeReturnSimInstances[simIndex].inventoryState[:, t], capacity, NB_PRODUCTS, addMmfeMean, addMmfeCovarianceMatrix, useFreeReturn=true)
            end
            extendPlaProdOut, extendProdRecourse, extendPlaSolveTime = extended_mmfe_with_pla(t, forecastStored[simIndex][:, :, t], extendMmfeSimInstances[simIndex].inventoryState[:, t], forecastEvolType, trueForecastEvolDis, capacity, T_BREAK, NB_PRODUCTS)
            extendFreeReturnPlaProdOut, extendFreeReturnProdRecourse, extendFreeReturnPlaSolveTime = extended_mmfe_with_pla(t, forecastStored[simIndex][:, :, t], extendFreeReturnSimInstances[simIndex].inventoryState[:, t], forecastEvolType, trueForecastEvolDis, capacity, T_BREAK, NB_PRODUCTS, useFreeReturn=true)

            ## Update simulation instances
            # MMFE - No returns
            update_instance_state!(t, mmfeProdOut, mmfePlaSolveTime, demandRealStored[simIndex][:, t], mmfeSimInstances[simIndex], NB_PRODUCTS)
            update_instance_state!(t, extendPlaProdOut, extendProdRecourse, extendPlaSolveTime, demandRealStored[simIndex][:, t], extendMmfeSimInstances[simIndex], NB_PRODUCTS)
            # MMFE - Free returns
            update_instance_state!(t, freeReturnMmfeProdOut, freeReturnMmfePlaSolveTime, demandRealStored[simIndex][:, t], freeReturnSimInstances[simIndex], NB_PRODUCTS)
            update_instance_state!(t, extendFreeReturnPlaProdOut, extendFreeReturnProdRecourse, extendFreeReturnPlaSolveTime, demandRealStored[simIndex][:, t], extendFreeReturnSimInstances[simIndex], NB_PRODUCTS)
            # Mis-specified MMFE model
            update_instance_state!(t, misspecMmfeProdOut, misspecMmfePlaSolveTime, demandRealStored[simIndex][:, t], misspecSimInstances[simIndex], NB_PRODUCTS)
            update_instance_state!(t, misspecFreeReturnMmfeProdOut, misspecFreeReturnMmfePlaSolveTime, demandRealStored[simIndex][:, t], misspecFreeReturnSimInstances[simIndex], NB_PRODUCTS)
        end
    end
    ## - KPI calculation -
    ## Calculate realised cost over simulation
    # MMFE - No returns
    mmfePlaCosts, mmfeSl, mmfeAggSl = kpis(mmfeSimInstances, backlogCost, demandRealStored, NB_PRODUCTS, NB_SIMULATIONS)
    extendMmfePlaCosts, extendMmfeSl, extendMmfeAggSl = kpis(extendMmfeSimInstances, backlogCost, demandRealStored, NB_PRODUCTS, NB_SIMULATIONS)
    # MMFE - Free returns
    freeReturnMmfePlaCosts, freeReturnMmfeSl, freeReturnMmfeAggSl = kpis(freeReturnSimInstances, backlogCost, demandRealStored, NB_PRODUCTS, NB_SIMULATIONS)
    extendFreeReturnPlaCosts, extendFreeReturnSl, extendFreeReturnAggSl = kpis(extendFreeReturnSimInstances, backlogCost, demandRealStored, NB_PRODUCTS, NB_SIMULATIONS)
    # Mis-specified MMFE model - Free returns but no recourse
    misspecMmfePlaCosts, misspecMmfeSl, misspecMmfeAggSl = kpis(misspecSimInstances, backlogCost, demandRealStored, NB_PRODUCTS, NB_SIMULATIONS)
    misspecFreeReturnMmfePlaCosts, misspecFreeReturnMmfeSl, misspecFreeReturnMmfeAggSl = kpis(misspecFreeReturnSimInstances, backlogCost, demandRealStored, NB_PRODUCTS, NB_SIMULATIONS)

    ## Store result in resultDataframe
    for s in 1:NB_SIMULATIONS
        push!(resultDataframe, [demandPattern forecastEvolType "mmfe" s mmfePlaCosts[s, 4] mmfeAggSl[s]])
        push!(resultDataframe, [demandPattern forecastEvolType "mmfe-FR" s freeReturnMmfePlaCosts[s, 4] freeReturnMmfeAggSl[s]])
        push!(resultDataframe, [demandPattern forecastEvolType "extended" s extendMmfePlaCosts[s, 4] extendMmfeAggSl[s]])
        push!(resultDataframe, [demandPattern forecastEvolType "extended-FR" s extendFreeReturnPlaCosts[s, 4] extendFreeReturnAggSl[s]])
        push!(resultDataframe, [demandPattern forecastEvolType "misspec" s misspecMmfePlaCosts[s, 4] misspecMmfeAggSl[s]])
        push!(resultDataframe, [demandPattern forecastEvolType "misspec-FR" s misspecFreeReturnMmfePlaCosts[s, 4] misspecFreeReturnMmfeAggSl[s]])
    end
end
elapsed = time() - start
println("\n---------- End of free return analysis, elapsed: $elapsed s. ----------")
# Export dataframe result to csv
CSV.write("output\\resultDataframe_free_return_cap_300.csv", resultDataframe)

## --- Plot results ---
# resultDataframe = CSV.read("output\\resultDataframe_free_return.csv")
addResultDf = filter(row -> row[:forecastEvolType] == "Additive", resultDataframe)
multResultDf = filter(row -> row[:forecastEvolType] == "Multiplicative", resultDataframe)

# -- Value of recourse --
# Read value of recourse from resultDataframe
addValOfRecourseMmfe, addValOfRecourseFreeReturns = value_of_recourse(addResultDf)
multValOfRecourseMmfe, multValOfRecourseFreeReturns = value_of_recourse(multResultDf)
boxplot([addValOfRecourseMmfe addValOfRecourseFreeReturns multValOfRecourseMmfe multValOfRecourseFreeReturns], title="Value of recourse", label = ["add-base" "add-fr" "mult-base" "mult-fr"])

# Export results to LaTeX table
resultTable = Matrix(undef, 2, 5)
resultTable[:, 1] = ["Value of recourse", "Difference"]
# Fill result table
resultTable[:, 2], resultTable[:, 3] = fill_column_recourse_table(addValOfRecourseMmfe, addValOfRecourseFreeReturns)
resultTable[:, 4], resultTable[:, 5] = fill_column_recourse_table(multValOfRecourseMmfe, multValOfRecourseFreeReturns)
latex_tabular(string("output\\FR-ValOfRecourse.tex"),
               Tabular("ccccc"),
               [Rule(:top),
               [" ","\\multicolumn{2}{c}{Additive MMFE}","\\multicolumn{2}{c}{Multiplicative MMFE}"],
               [" \\cmidrule(lr){2-3} \\cmidrule(lr){4-5} ","No returns", "Free returns", "No returns", "Free returns"],
               Rule(:mid),
               resultTable,
               Rule(:bottom)])

# Value of free returns
addAvgCostsNoReturns, addAvgCostsFR, addValFreeReturnsPLA, addValFreeReturnsExt, addValFreeReturnsMisSpec = value_of_free_returns(addResultDf)
multAvgCostsNoReturns, multAvgCostsFR, multValFreeReturnsPLA, multValFreeReturnsExt, multValFreeReturnsMisSpec = value_of_free_returns(multResultDf)

# Export results to LaTeX table
resultTable = Matrix(undef, 3, 7)
resultTable[1:3, 1] = ["Avg. costs (no returns) \$\\bar{z}^{(\\text{nr})}\$",
                       "Avg. costs (free returns) \$\\bar{z}^{(\\text{fr})}\$",
                       "Value of returns \$\\bar{\\nu}\$"]
resultTable[1, 2:4] = addAvgCostsNoReturns
resultTable[1, 5:7] = multAvgCostsNoReturns
resultTable[2, 2:4] = addAvgCostsFR
resultTable[2, 5:7] = multAvgCostsFR
# Fill result table
resultTable[3, 2:4] = format_val_of_returns(addValFreeReturnsPLA, addValFreeReturnsExt, addValFreeReturnsMisSpec)
resultTable[3, 5:7] = format_val_of_returns(multValFreeReturnsPLA, multValFreeReturnsExt, multValFreeReturnsMisSpec)
latex_tabular(string("output\\FR-ValOfReturns.tex"),
               Tabular("lcccccc"),
               [Rule(:top),
               [" ","\\multicolumn{3}{c}{Additive MMFE}","\\multicolumn{3}{c}{Multiplicative MMFE}"],
               [" \\cmidrule(lr){2-4} \\cmidrule(lr){5-7} ","Add. PLA", "Add. Ext.", "Mult. PLA", "Mult. PLA", "Mult. Ext", "Add. PLA"],
               ["","(correct)", "(correct)", "(mismatched)","(correct)", "(correct)", "(mismatched)"],
               Rule(:mid),
               resultTable,
               Rule(:bottom)])
