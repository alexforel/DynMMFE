function latin_hypercube_strata(nbSamples)
    ## latin_hypercube_strata(nbSamples)
    # Provide nbSamples uniformly distributed over the N-dimensions space. Samples are repetitively created and samples that have the minimum distance to all other samples are deleted. This allows to cover the N-dimensions space more uniformly than straight sampling.
    # Input:
    #   - nbSamples: scalar, number of scenarios
    # Output:
    #   - unifSamples: Matrix[nbDimensions ; nbSamples], uniform samples in the N dimensions

    println("Generating Latin Hypercube strata for $nbSamples samples.")
    # Step1: generate N x ML uniform random numbers
    N = NB_PRODUCTS * T # number of dimensions
    if nbSamples < 5
        M = 20
    elseif nbSamples < 10
        M = 10
    else
        M = 5
    end
    L = nbSamples
    unifSamples = rand(N, M * L)
    # Step2: calculate euclidian distance
    dimension = M * L
    while dimension > L
        euclidianDistance = zeros(dimension, dimension)
        smallestAvergaeDistance = zeros(dimension)
        for i in 1:dimension
            for j in 1:dimension
                if i == j
                    euclidianDistance[i, j] = Inf
                else
                    euclidianDistance[i, j] = sqrt(sum(unifSamples[:, i] - unifSamples[:, j]).^2)
                end
            end
            # Find min distance
            minDist, minIndex = findmin(euclidianDistance[i, :])
            secondMinDist = findmin(euclidianDistance[i, 1:end .!= minIndex])[1]
            smallestAvergaeDistance[i] = mean([minDist secondMinDist])
        end
        # Delete sample with min distance and iterate counter
        minSample = findmin(smallestAvergaeDistance)[2]
        unifSamples = unifSamples[:, 1:end .!= minSample]
        dimension = dimension - 1
    end

    return unifSamples
end

function latin_hypercube_sampling_with_multidimensional_uniformity(latinStrata, inputMultivDistribution, nbSamples)
    ## latin_hypercube_sampling_with_multidimensional_uniformity(latinStrata, inputMultivDistribution, nbSamples)
    # Implementation of multivariate correlated Latin Hypercube. Adapted from Deutsch, J. L., & Deutsch, C. V. (2012). Latin hypercube sampling with multidimensional uniformity. Journal of Statistical Planning and Inference, 142(3), 763-772.
    # Input:
    #   - latinStrata: Matrix[N ; nbSamples], uniform samples
    #   - inputMultivDistribution: input multivariate distribution to be sampled
    #   - nbSamples
    # Output:
    #   - correlatedSamples

    # Read input data
    covarianceMatrix = cov(inputMultivDistribution)
    correlationMatrix = cor(inputMultivDistribution)
    # Define parameters and initialise matrices
    N = NB_PRODUCTS * T
    latinUniformSamples = zeros(N, nbSamples)
    latinSamples = zeros(N, nbSamples)

    # Following steps to correlate samples have been replaced by Cholesky multiplication at the end of procedure
    # # Step 6a: transform N x L values to Gaussian units using the standard Gaussian cdf
    # unitGaussians = Distributions.quantile.(Normal(), latinStrata)
    # # Step 6b:multiply each realisation i in L by the L matrix from holesky decomposition
    # correlatedGaussians = cholesky(correlationMatrix).L * unitGaussians
    # # Step 6c: backtransform to uniform distributions using CDF
    # correlatedStrata = Distributions.cdf.(Normal(), correlatedGaussians)

    for n in 1:N
        tempDis = Normal(mean(inputMultivDistribution)[n], sqrt(covarianceMatrix[n,n]))
        # Step 7: rank inputs
        sortedIndicesOfStrata = sortperm(latinStrata[n,:])
        for l in 1:nbSamples
            # Step 8: generate uniform numbers for the L strata
            latinUniformSamples[n, l] = rand(Uniform((sortedIndicesOfStrata[l] - 1)/nbSamples, sortedIndicesOfStrata[l]/nbSamples))
            # Step 9: sample the cumulative distribution of n with uniform numbers from Step 8
            latinSamples[n, l] = Distributions.quantile(tempDis, latinUniformSamples[n, l])
        end
    end

    # Correlate samples with cholesky L matrix
    correlatedSamples = cholesky(correlationMatrix).L * latinSamples

    return correlatedSamples
end
