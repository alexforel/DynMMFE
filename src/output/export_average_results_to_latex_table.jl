function export_average_results_to_latex_table(addAvgCost, addRelCost, addStatSignMatrix, multAvgCost, multRelCost, multStatSignMatrix)
    ## export_average_results_to_latex_table(inputAverageCostMatrix, inputRelativeCostMatrix)
    # Create a table providing abolute and relative costs of all simulation models compared to deterministic cost. All siulation parameters are presented from left to right in columns.
    # Input:
    #   - addAvgCost, addRelCost, addStatSignMatrix, multAvgCost, multRelCost, multStatSignMatrix
    # Output:
    #   - nothing

    # Format all results in relevant matrix
    addResultTable = format_results_in_table(addAvgCost, addRelCost, addStatSignMatrix, "Additive")
    multResultTable = format_results_in_table(multAvgCost, multRelCost, multStatSignMatrix, "Multiplicative")
    # Derive average performances over all instances and format it
    addLastRow = last_average_row_to_stor(addAvgCost, addStatSignMatrix)
    multLastRow = last_average_row_to_stor(multAvgCost, multStatSignMatrix)
    # Create table and save it as .tex file
    latex_tabular(string("output\\addAverageResultsSynthDataTable.tex"),
                   Tabular("cccccccc"),
                   [Rule(:top),
                   ["Demand", "Uncertainty", "\$cap\$", "Det.", "DD-stochastic", "Additive PLA", "Multiplicative PLA", "Ext. Add. PLA"],
                   Rule(:mid),
                   addResultTable,
                   Rule(:mid),
                   addLastRow,
                   Rule(:bottom)])
   latex_tabular(string("output\\multAverageResultsSynthDataTable.tex"),
                  Tabular("cccccccc"),
                  [Rule(:top),
                  ["Demand", "Uncertainty", "\$cap\$", "Det.", "DD-stochastic", "Additive PLA", "Multiplicative PLA", "Ext. Mult. PLA"],
                  Rule(:mid),
                  multResultTable,
                  Rule(:mid),
                  multLastRow,
                  Rule(:bottom)])
    return nothing
end

function format_results_in_table(inputAverageCostMatrix, inputRelativeCostMatrix, statSignMatrix, mmfeInputString)
    # Define matrix to store
    formattedResultTable = Matrix(undef, size(inputAverageCostMatrix)[1], 3 + 5)
    for demP in demandPatternVector
        for unCo in uncertaintyCoeffVector
            for cap in capacityVector
                # Find row at which to store results
                a = findfirst(isequal(demP), demandPatternVector)
                b = findfirst(isequal(unCo), uncertaintyCoeffVector)
                c = findfirst(isequal(cap), capacityVector)
                rowToStore = length(capacityVector) * length(uncertaintyCoeffVector) * (a-1) + length(capacityVector) * (b-1) + c
                ## Format left side of table matrix
                formattedResultTable[rowToStore, 1] = ""
                if (rowToStore % (length(uncertaintyCoeffVector) * length(capacityVector)) == 1)
                    formattedResultTable[rowToStore, 1] =  string("\\multirow{4}{*}{",demP,"}")
                else
                    formattedResultTable[rowToStore, 1] = ""
                end
                if (rowToStore % length(capacityVector) == 1)
                    if unCo == 1
                        formattedResultTable[rowToStore, 2] =  string("\\multirow{2}{*}{Low}")
                    elseif unCo == 2
                        formattedResultTable[rowToStore, 2] =  string("\\multirow{2}{*}{Medium}")
                    elseif unCo == 3
                        formattedResultTable[rowToStore, 2] =  string("\\multirow{2}{*}{High}")
                    end
                else
                    formattedResultTable[rowToStore, 2] = ""
                end
                # Store numerical results
                floatVector = round.([inputAverageCostMatrix[rowToStore, 1] inputAverageCostMatrix[rowToStore, 2] inputRelativeCostMatrix[rowToStore, 2] inputAverageCostMatrix[rowToStore, 3] inputRelativeCostMatrix[rowToStore, 3] inputAverageCostMatrix[rowToStore, 4] inputRelativeCostMatrix[rowToStore, 4] inputAverageCostMatrix[rowToStore, 5] inputRelativeCostMatrix[rowToStore, 5]], digits = 1)
                formattedResultTable[rowToStore, 3:4] = [cap floatVector[1]]
                # Parse results with statistical significance
                if statSignMatrix[rowToStore, 1] == 1.0
                    formattedResultTable[rowToStore, 5] = string(floatVector[2], "\\, (", floatVector[3], "\\% (*))")
                else
                    formattedResultTable[rowToStore, 5] = string(floatVector[2], "\\, (", floatVector[3], "\\%)")
                end
                if statSignMatrix[rowToStore, 2] == 1.0
                    formattedResultTable[rowToStore, 6] = string(floatVector[4], "\\, (", floatVector[5], "\\% (*))")
                else
                    formattedResultTable[rowToStore, 6] = string(floatVector[4], "\\, (", floatVector[5], "\\%)")
                end
                if statSignMatrix[rowToStore, 3] == 1.0
                    formattedResultTable[rowToStore, 7] = string(floatVector[6], "\\, (", floatVector[7], "\\% (*))")
                else
                    formattedResultTable[rowToStore, 7] = string(floatVector[6], "\\, (", floatVector[7], "\\%)")
                end
                if statSignMatrix[rowToStore, 4] == 1.0
                    formattedResultTable[rowToStore, 8] = string(floatVector[8], "\\, (", floatVector[9], "\\% (*))")
                else
                    formattedResultTable[rowToStore, 8] = string(floatVector[8], "\\, (", floatVector[9], "\\%)")
                end

            end
        end
    end
    return formattedResultTable
end

function last_average_row_to_stor(inputAverageCostMatrix, statSignMatrix)
    ## Define last row which is average over all simulation settings

    # Calculate all costs to store
    detCost = round(mean(inputAverageCostMatrix[:, 1]), digits = 1)
    ddstochCost = round(mean(inputAverageCostMatrix[:, 2]), digits = 1)
    ddstochRelCost = round(ddstochCost/detCost*100, digits = 1)
    addCost = round(mean(inputAverageCostMatrix[:, 3]), digits = 1)
    addRelCost = round(addCost/detCost*100, digits = 1)
    multCost = round(mean(inputAverageCostMatrix[:, 4]), digits = 1)
    multRelCost = round(multCost/detCost*100, digits = 1)
    extCost = round(mean(inputAverageCostMatrix[:, 5]), digits = 1)
    extRelCost = round(extCost/detCost*100, digits = 1)

    # Store costs in vector
    lastRowToStore = Matrix(undef, 1, 1 + 5)
    lastRowToStore[1] = MultiColumn(3, :c, "Average")
    lastRowToStore[2] = detCost
    if statSignMatrix[end, 1] == 1.0
        lastRowToStore[3] = string(ddstochCost, "\\, (", ddstochRelCost, "\\% (*))")
    else
        lastRowToStore[3] = string(ddstochCost, "\\, (", ddstochRelCost, "\\%)")
    end
    if statSignMatrix[end, 2] == 1.0
        lastRowToStore[4] = string(addCost, "\\, (", addRelCost, "\\% (*))")
    else
        lastRowToStore[4] = string(addCost, "\\, (", addRelCost, "\\%)")
    end
    if statSignMatrix[end, 3] == 1.0
        lastRowToStore[5] = string(multCost, "\\, (", multRelCost, "\\% (*))")
    else
        lastRowToStore[5] = string(multCost, "\\, (", multRelCost, "\\%)")
    end
    if statSignMatrix[end, 4] == 1.0
        lastRowToStore[6] = string(extCost, "\\, (", extRelCost, "\\% (*))")
    else
        lastRowToStore[6] = string(extCost, "\\, (", extRelCost, "\\%)")
    end

    return lastRowToStore
end
