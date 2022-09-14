# -- Illustrate results on MMFE and cumulative demand. --
# Reproduce figures from paper:
# - Figure 2. Evolution of variance with correlation coefficient for the (a) additive and (b) multiplicative MMFE.
# - Figure 3. PLA of expected inventory and backlog for demand following (a) normal distribution and (b)log-normal distribution
# - Figure 7. Mean demand for (a) stationary, (b) random, and (c) seasonal patterns over simulation of 8 periods

## Initialisation
# Set working directory
cd("C:/Users/Alexandre//Documents/Code/DynMMFE")

using Random, Distributions, Statistics
using QuadGK
using DataFrames, CSV

# Load functions
include(raw"src\analysis\exact_and_pla_folf.jl")
include(raw"src\analysis\cumulated_demand_variance_correlation_analysis.jl")
include(raw"src\simpson_integration.jl")
include(raw"src\demand\demand_distribution_pattern.jl")
include(raw"src\mmfe\fenton_wilkinson_approx_sum_log_normal.jl")
include(raw"src\mmfe\cumulative_demand_distribution_matrix_from_additive_mmfe.jl")
include(raw"src\mmfe\cumulative_demand_distribution_matrix_from_multiplicative_mmfe.jl")
include(raw"src\folf\first_order_loss_function.jl")
include(raw"src\folf\folf_linearisation.jl")
# Types
include(raw"src\types\PiecewiseLinear.jl")

# Piecewise-linear approximation of the inverse of the standard folf
xLim = 5
first_order_loss_function_standard_normal(x) = quadgk(t -> (1 - cdf(Normal(), t)), x, Inf)[1]
epsilonRange = collect(range(-xLim, stop = xLim, length = 2000))
pwlaStandardFolf = PiecewiseLinear(epsilonRange, first_order_loss_function_standard_normal)

## 1/ ----  Evolution of variance with correlation coefficient
# For one product, over two periods: CD_{k=1, t=2}
nbProducts = 1 ;T = 2
nbCor = 11
correlationVector = collect(range(-0.999, stop = 0.999, length =nbCor))
t1 = 1 ; t2 = 2
# Additive MMFE
sigma1 = 20 ; sigma2 = 20
# Multiplicative model
nbSamples = 100000
multvar = 0.040
# Variance of cumulative demand
cumDemandVarianceAdd50, cumDemandVarianceFw50, cumDemandVarianceSamples50 = cumulated_demand_variance_correlation_analysis(50)
cumDemandVarianceAdd100, cumDemandVarianceFw100, cumDemandVarianceSamples100 = cumulated_demand_variance_correlation_analysis(100)
cumDemandVarianceAdd150, cumDemandVarianceFw150, cumDemandVarianceSamples150 = cumulated_demand_variance_correlation_analysis(150)
# Export results to csv file
varianceEvolutionDataFrame = DataFrame([correlationVector cumDemandVarianceAdd50 cumDemandVarianceFw50 cumDemandVarianceFw100 cumDemandVarianceFw150], :auto)
CSV.write("output\\csv\\varianceEvolution.csv", varianceEvolutionDataFrame)

## 2/ ---- PLA of expected inventory and backlog ---
T = 1 ; VERBOSE = 0 ; nbSegments = 6
meanDemand = 100
# NORMAL DISTRIBUTION
addVar = 400
normalDistribution = Normal(meanDemand, sqrt(addVar))
xVect, normFolf, normComplfFolf, uRange, normlinFolf, normlinComplFolf = exact_and_pla_folf(normalDistribution)
# LOGNORMAL DISTRIBUTION
multVar = 0.04
logNormalDistribution = LogNormal(log(meanDemand)-0.5*multVar, sqrt(multVar))
xVect, lognormFolf, lognormComplfFolf, uRange, lognormlinFolf, lognormlinComplFolf = exact_and_pla_folf(logNormalDistribution)
# Export to csv tables
folfDf = DataFrame([xVect normFolf normComplfFolf lognormFolf lognormComplfFolf], :auto)
CSV.write("output\\csv\\folfDf.csv", folfDf)
linFolfDf = DataFrame([uRange[1,1,:] normlinFolf normlinComplFolf lognormlinFolf lognormlinComplFolf], :auto)
CSV.write("output\\csv\\linFolfDf.csv", linFolfDf)

## 3/ ----- Mean demand for (a) stationary, (b) random, and (c) seasonal patterns -----
T = 6 # prediction horizon, and period of seasonal demand
# Visualise demand patterns
staticDistributions = demand_distribution_pattern("Stationary", 12)
randDistributions = demand_distribution_pattern("Random", 12)
seasonalDistributions = demand_distribution_pattern("Seasonal", 12)
# Export mean demands to csv for plot in latex
patternDf = DataFrame([1:12 mean.(staticDistributions) rand.(randDistributions) mean.(seasonalDistributions)], :auto)
CSV.write("output\\demandPatterns.csv", patternDf)
