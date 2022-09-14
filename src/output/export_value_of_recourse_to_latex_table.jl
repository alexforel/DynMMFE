function export_value_of_recourse_to_latex_table(inputDf, inputString, fileNameString)
    ## export_value_of_recourse_to_latex_table(inputDf, inputString, fileNameString)
    # Create a table providing abolute and relative costs of all simulation models compared to deterministic cost. All siulation parameters are presented from left to right in columns.
    # Input:
    #   - inputDf, inputString, fileNameString
    # Output:
    #   - nothing

    firstRow = ["\\multirow{2}{*}{Demand}"; "\\multirow{2}{*}{Uncertainty}"; "\\multirow{2}{*}{Capacity}"; "\\multirow{2}{*}{PLA cost}";"\\multicolumn{5}{c}{Extended PLA with recourse}"]
    secondRow = ["";"";"";""; "Absolute cost"; "Relative avg."; "First quart."; "Median"; "Third quart."]
    formattedResults = Matrix(undef, size(addAvgCost)[1], 3 + 6)
    for demP in demandPatternVector
        tempDf = filter(row -> row[:demandPattern] == demP, inputDf)
        for unCo in uncertaintyCoeffVector
            tempDf2 = filter(row -> row[:uncertaintyCoeff] == unCo, tempDf)
            for cap in capacityVector
                # Filter dataframe to have only cost of no recourse and extended
                tempDf3 = filter(row -> row[:capacity] == cap, tempDf2)
                noRecourseCost = filter(row -> row[:model] == inputString, tempDf3)[:cost]
                extendedCost = filter(row -> row[:model] == "extended", tempDf3)[:cost]
                # Fid row at which to store results
                a = findfirst(isequal(demP), demandPatternVector)
                b = findfirst(isequal(unCo), uncertaintyCoeffVector)
                c = findfirst(isequal(cap), capacityVector)
                rowToStore = length(capacityVector) * length(uncertaintyCoeffVector) * (a-1) + length(capacityVector) * (b-1) + c
                # Store relevant values
                if (rowToStore % (length(uncertaintyCoeffVector) * length(capacityVector)) == 1)
                    formattedResults[rowToStore, 1] =  string("\\multirow{4}{*}{",demP,"}")
                else
                    formattedResults[rowToStore, 1] = ""
                end
                if (rowToStore % length(capacityVector) == 1)
                    if unCo == 1
                        formattedResults[rowToStore, 2] =  string("\\multirow{2}{*}{Low}")
                    elseif unCo == 2
                        formattedResults[rowToStore, 2] =  string("\\multirow{2}{*}{Medium}")
                    elseif unCo == 3
                        formattedResults[rowToStore, 2] =  string("\\multirow{2}{*}{High}")
                    end
                else
                    formattedResults[rowToStore, 2] = ""
                end
                formattedResults[rowToStore,3] = cap
                formattedResults[rowToStore,4] = round(mean(noRecourseCost), digits = 1)
                formattedResults[rowToStore,5] = round(mean(extendedCost), digits = 1)
                ## Relative cost
                relativeCost = extendedCost ./ noRecourseCost * 100
                statSign = pvalue(OneSampleTTest(relativeCost, 100)) <= 0.05
                # Store avg relative cost with statistical significance
                avgRelCost = round(mean(relativeCost), digits = 1)
                if statSign
                    formattedResults[rowToStore,6] = string(avgRelCost, "\\% (*)")
                else
                    formattedResults[rowToStore,6] = string(avgRelCost, "\\%")
                end
                firstQuart = round(quantile(relativeCost, 0.25), digits = 1)
                formattedResults[rowToStore,7] = string(firstQuart, "\\%")
                medianC = round(median(relativeCost), digits = 1)
                formattedResults[rowToStore,8] = string(medianC, "\\%")
                thirdQuart = round(quantile(relativeCost, 0.75), digits = 1)
                formattedResults[rowToStore,9] = string(thirdQuart, "\\%")
            end
        end
    end
    # Define last row which is average over all simulation settings
    noRecourseCost = filter(row -> row[:model] == inputString, inputDf)[:cost]
    extendedCost = filter(row -> row[:model] == "extended", inputDf)[:cost]
    relativeCost = 100 * (extendedCost./noRecourseCost)
    # Assess statistical significance of average
    statSign = pvalue(OneSampleTTest(relativeCost, 100)) <= 0.05

    lastRowToStore = Matrix(undef, 1, 1 + 6)
    lastRowToStore[1] = MultiColumn(3, :c, "All settings")
    lastRowToStore[2] = round(mean(noRecourseCost), digits = 1)
    lastRowToStore[3] = round(mean(extendedCost), digits = 1)
    # Store avg relative cost with statistical significance
    avgRelCost = round(mean(relativeCost), digits = 1)
    if statSign
        lastRowToStore[4] = string(avgRelCost, "\\% (*)")
    else
        lastRowToStore[4] = string(avgRelCost, "\\%")
    end
    firstQuart = round(quantile(relativeCost, 0.25), digits = 1)
    lastRowToStore[5] = string(firstQuart, "\\%")
    medianC = round(median(relativeCost), digits = 1)
    lastRowToStore[6] = string(medianC, "\\%")
    thirdQuart = round(quantile(relativeCost, 0.75), digits = 1)
    lastRowToStore[7] = string(thirdQuart, "\\%")

    latex_tabular(string("output\\",fileNameString,".tex"), Tabular("ccccccccc"),
                   [Rule(:top), firstRow, secondRow, Rule(:mid), formattedResults, Rule(), lastRowToStore, Rule(:bottom)])

    return nothing
end
