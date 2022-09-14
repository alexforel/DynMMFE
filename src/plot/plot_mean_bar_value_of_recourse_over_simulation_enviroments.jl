function plot_mean_bar_value_of_recourse_over_simulation_enviroments(inputResultDf, noRecourseModelString)
    ## plot_mean_bar_value_of_recourse_over_simulation_enviroments(inputResultDf, noRecourseModelString)

    plotList = []
    for demP in demandPatternVector
        tempDf = filter(row -> row[:demandPattern] == demP, inputResultDf)
        for c in capacityVector
            tempDf2 = filter(row -> row[:capacity] == c, tempDf)
            for u in uncertaintyCoeffVector
                tempDf3 = filter(row -> row[:uncertaintyCoeff] == u, tempDf2)
                noRecourseCost = filter(row -> row[:model] == noRecourseModelString, tempDf3)[:cost]
                extCost = filter(row -> row[:model] == "extended", tempDf3)[:cost]
                # Bar plot of mean value
                push!(plotList, bar([mean(extCost./noRecourseCost)], legend = false, title = "$demP, cap = $c, u = $u"))
            end
        end
    end
    # Output
    outPlot = plot(plotList..., size = (1200, 800), fontfamily = "Times", xticks = (1:2, forecastEvolTypeVector), ylims = (0.9, 1.1))
    return outPlot
end
