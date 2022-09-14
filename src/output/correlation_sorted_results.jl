function correlation_sorted_results(inputDf)
    ## correlation_sorted_results(inputDf)
    # Filter result data frame to the relevant sensitivity instance. Calculate confidence interval on absolute and relative costs. Perform hypothesis tests.
    # Input:
    #   - inputDf
    # Output:
    #   - mmfeAvgCost, mmfeConf, extMmfeAvgCost, extendConf, relativeCost, statSign

    mmfeAvgCost = zeros(length(prodCorrVector), length(timeCorrVector))
    extMmfeAvgCost = zeros(length(prodCorrVector), length(timeCorrVector))
    mmfeConf = zeros(length(prodCorrVector), length(timeCorrVector))
    extendConf = zeros(length(prodCorrVector), length(timeCorrVector))
    relativeCost = zeros(length(prodCorrVector), length(timeCorrVector))
    statSign = zeros(length(prodCorrVector), length(timeCorrVector))
    for pc in 1:length(prodCorrVector)
        tempDf = filter(row -> row[:prodCorr] == prodCorrVector[pc], inputDf)
        for tc in 1:length(timeCorrVector)
            tempDf2 = filter(row -> row[:timeCorr] == timeCorrVector[tc], tempDf)
            # Filter df, get results and confidence interval
            tempMmfeDf = filter(row -> row[:model] == "mmfe", tempDf2)
            mmfeSamples = tempMmfeDf[!, :cost]
            mmfeAvgCost[pc, tc] = mean(mmfeSamples)
            mmfeConf[pc, tc] = 1.96 * std(mmfeSamples) / sqrt(length(mmfeSamples))
            tempExtMmfeDf = filter(row -> row[:model] == "extended", tempDf2)
            extSamples = tempExtMmfeDf[!, :cost]
            extMmfeAvgCost[pc, tc] = mean(extSamples)
            extendConf[pc, tc] = 1.96 * std(extSamples) / sqrt(length(extSamples))
            ## Relative cost and its statistical significance
            relativeCostSamples = extSamples ./ mmfeSamples * 100
            relativeCost[pc, tc] = mean(relativeCostSamples)
            # Statistical significance of cost difference
            statSign[pc, tc] = pvalue(OneSampleTTest(relativeCostSamples, 100)) <= 0.05
        end
    end

    return  mmfeAvgCost, mmfeConf, extMmfeAvgCost, extendConf, relativeCost, statSign
end
