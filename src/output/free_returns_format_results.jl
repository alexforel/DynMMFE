function value_of_recourse(inputDf)
    # Read value of recourse from resultDataframe
    mmfeCost = filter(row -> row[:model] == "mmfe", inputDf)[!, :cost]
    mmfeFRCost = filter(row -> row[:model] == "mmfe-FR", inputDf)[!, :cost]
    extMmfeCost = filter(row -> row[:model] == "extended", inputDf)[!, :cost]
    extMmfeFRCost = filter(row -> row[:model] == "extended-FR", inputDf)[!, :cost]

    valOfRecourseMmfe = extMmfeCost./mmfeCost * 100
    valOfRecourseFreeReturns = extMmfeFRCost./mmfeFRCost * 100

    return valOfRecourseMmfe, valOfRecourseFreeReturns
end

function column_percentage_point_diff(valMmfe, valMethod2)
    valDiff = valMethod2 - valMmfe
    statSign = pvalue(OneSampleTTest(valDiff, 0.0)) <= 0.05
    meanValDiff = round(mean(valDiff), digits = 2)
    meanValMethod2 = round(mean(valMethod2), digits = 1)
    if statSign
        col2 = [string(meanValMethod2, "\\%"),
                             string(meanValDiff, " pp (*)")]
    else
        col2 = [string(meanValMethod2, "\\%"),
                             string(meanValDiff, " pp")]
    end

    return col2
end

function fill_column_recourse_table(valOfRecourseMmfe, valOfRecourseFreeReturns)
    meanValOfRecourseMmfe = round(mean(valOfRecourseMmfe), digits = 1)
    col1 = [string(meanValOfRecourseMmfe, "\\%"), ""]
    # Difference in Value of Recourse: statistically significant?
    col2 = column_percentage_point_diff(valOfRecourseMmfe, valOfRecourseFreeReturns)
    return col1, col2
end

function value_of_free_returns(inputDf)
    # Read value of recourse from resultDataframe
    mmfeCost = filter(row -> row[:model] == "mmfe", inputDf)[!, :cost]
    mmfeFRCost = filter(row -> row[:model] == "mmfe-FR", inputDf)[!, :cost]
    extMmfeCost = filter(row -> row[:model] == "extended", inputDf)[!, :cost]
    extMmfeFRCost = filter(row -> row[:model] == "extended-FR", inputDf)[!, :cost]
    misspecMmfeCost = filter(row -> row[:model] == "misspec", inputDf)[!, :cost]
    misspecMmfeFRCost = filter(row -> row[:model] == "misspec-FR", inputDf)[!, :cost]
    # Returns average sim costs
    avgPlaNoReturns = round(mean(mmfeCost), digits = 1)
    avgPlaFR = round(mean(mmfeFRCost), digits = 1)
    avgExtNoReturns = round(mean(extMmfeCost), digits = 1)
    avgExtFR = round(mean(extMmfeFRCost), digits = 1)
    avgMispecNoReturns = round(mean(misspecMmfeCost), digits = 1)
    avgMispecFR = round(mean(misspecMmfeFRCost), digits = 1)
    avgCostsNoReturns = [avgPlaNoReturns, avgExtNoReturns, avgMispecNoReturns]
    avgCostsFR = [avgPlaFR, avgExtFR, avgMispecFR]
    # Calculate value of free returns
    valFreeReturnsPLA = mmfeFRCost./mmfeCost * 100
    valFreeReturnsExt = extMmfeFRCost./extMmfeCost * 100
    valFreeReturnsMisSpec = misspecMmfeFRCost./misspecMmfeCost * 100

    return avgCostsNoReturns, avgCostsFR, valFreeReturnsPLA, valFreeReturnsExt, valFreeReturnsMisSpec
end


function format_val_of_returns(valFreeReturnsPLA, valFreeReturnsExt, valFreeReturnsMisSpec)
    # Calculate average values of free returns for all models
    meanValReturnsPla = round(mean(valFreeReturnsPLA), digits = 1)
    statSign = pvalue(OneSampleTTest(valFreeReturnsPLA, 100)) <= 0.05
    if statSign
        val1 = string(meanValReturnsPla, "\\% (*)")
    else
        val1 = string(meanValReturnsPla, "\\%")
    end

    meanValReturnsExt = round(mean(valFreeReturnsExt), digits = 1)
    statSign = pvalue(OneSampleTTest(valFreeReturnsExt, 100)) <= 0.05
    if statSign
        val2 = string(meanValReturnsExt, "\\% (*)")
    else
        val2 = string(meanValReturnsExt, "\\%")
    end

    meanValReturnsMisSpec = round(mean(valFreeReturnsMisSpec), digits = 1)
    statSign = pvalue(OneSampleTTest(valFreeReturnsMisSpec, 100)) <= 0.05
    if statSign
        val3 = string(meanValReturnsMisSpec, "\\% (*)")
    else
        val3 = string(meanValReturnsMisSpec, "\\%")
    end
    return [val1, val2, val3]
end
