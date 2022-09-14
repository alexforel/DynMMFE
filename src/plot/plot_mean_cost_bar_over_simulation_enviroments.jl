function plot_mean_cost_bar_over_simulation_enviroments(inputResultDf)
    ## plot_mean_cost_bar_over_simulation_enviroments(inputResultDf)
    # Plot avergae costs as a bar chart. Each parameter configuration is a subplot.
    # Input:
    #   - inputResultDf
    # Output:
    #   - outPlot

    plotList = []
    for demP in demandPatternVector
        tempDf = filter(row -> row[:demandPattern] == demP, inputResultDf)
        for c in capacityVector
            tempDf2 = filter(row -> row[:capacity] == c, tempDf)
            for u in uncertaintyCoeffVector
                tempDf3 = filter(row -> row[:uncertaintyCoeff] == u, tempDf2)
                detCost = filter(row -> row[:model] == "deterministic", tempDf3)[!, :cost]
                ddstochCost = filter(row -> row[:model] == "dd-stoch", tempDf3)[!, :cost]
                addCost = filter(row -> row[:model] == "additive", tempDf3)[!, :cost]
                multCost = filter(row -> row[:model] == "multiplicative", tempDf3)[!, :cost]
                # Mean bar
            	push!(plotList, bar([mean(ddstochCost./detCost),mean(addCost./detCost), mean(multCost./detCost)], legend = false, title = "$demP, cap = $c, u = $u"))
            end
        end
    end

    # Output
    outPlot = plot(plotList..., size = (1200, 800), fontfamily = "Times", xticks = (1:2, forecastEvolTypeVector))
    return outPlot
end
