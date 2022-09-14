function fenton_wilkinson_approx_sum_log_normal(logMeanVector, logSigmaCovar)
    ## fenton_wilkinson_approx_sum_log_normal(logMeanVector, logSigmaCovar)
    # Implemented for abitrary correlated log normal RVs as in Abu-Dayya, A. A.,
    # & Beaulieu, N. C. (1994). Outage probabilities in the presence of correlated
    # lognormal interferers. IEEE Transactions on Vehicular Technology, 43(1), 164-173.
    # Input:
    #   - logMeanVector
    #   - logSigmaCovar
    # Output:
    #   - approxLogNormalDistribution: Distribution, approximate log-normal (demand) distribution

    nbDistribution = length(logMeanVector)

    # Abu-Dayaa - Fenton-Wilkinson approximation
    u1 = sum(exp(logMeanVector[j] + 0.5 * logSigmaCovar[j,j]) for j in 1:nbDistribution)
    u2 = sum(exp(2*logMeanVector[j] + 2 * logSigmaCovar[j,j]) for j in 1:nbDistribution) + 2 * sum(sum(exp(logMeanVector[i] + logMeanVector[j]) * exp(0.5 * (logSigmaCovar[i,i] + logSigmaCovar[j,j] + 2 * logSigmaCovar[i,j]))  for j in (i+1):(nbDistribution)) for i in 1:(nbDistribution-1))

    # Determine two first log-moments of the sum of log normal RVs
    muZ = 2 * log(u1) - 0.5 * log(u2)
    sigmaZ = sqrt(log(u2) - 2 * log(u1))

    # Output
    approxLogNormalDistribution = LogNormal(muZ, sigmaZ)

    return approxLogNormalDistribution
end
