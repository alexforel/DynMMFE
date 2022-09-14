function plot_boxplot_value_of_recourse(addResultDf, multResultDf)
    ## plot_boxplot_value_of_recourse(inputResultDf)

    ## Plot 1 / Low capacity
    # Additive
    addNoRecourseCost = filter(row -> (row[:model] == "additive") & (row[:capacity] == 250), addResultDf)[:cost]
    addExtCost = filter(row ->(row[:model] == "extended") & (row[:capacity] == 250), addResultDf)[:cost]
    # Multiplicative
    multNoRecourseCost = filter(row -> (row[:model] == "multiplicative") & (row[:capacity] == 250), multResultDf)[:cost]
    multExtCost = filter(row -> (row[:model] == "extended") & (row[:capacity] == 250), multResultDf)[:cost]
    addAvgValue = mean(addExtCost./addNoRecourseCost) * 100
    multAvgValue = mean(multExtCost./multNoRecourseCost) * 100
    println("Low capacity: Average cost reduction of recourse is ADD: $addAvgValue , MULT: $multAvgValue.")
    plot1 = boxplot([100 * addExtCost ./ addNoRecourseCost 100 * multExtCost ./ multNoRecourseCost], outliers = false, title = "(a) Low capacity")

    ## Plot 2 / High capacity
    # Additive
    addNoRecourseCost = filter(row -> (row[:model] == "additive") & (row[:capacity] == 500), addResultDf)[:cost]
    addExtCost = filter(row ->(row[:model] == "extended") & (row[:capacity] == 500), addResultDf)[:cost]
    # Multiplicative
    multNoRecourseCost = filter(row -> (row[:model] == "multiplicative") & (row[:capacity] == 500), multResultDf)[:cost]
    multExtCost = filter(row -> (row[:model] == "extended") & (row[:capacity] == 500), multResultDf)[:cost]
    addAvgValue = mean(addExtCost./addNoRecourseCost) * 100
    multAvgValue = mean(multExtCost./multNoRecourseCost) * 100
    println("High capacity: Average cost reduction of recourse is ADD: $addAvgValue , MULT: $multAvgValue.")
    plot2 = boxplot([100 * addExtCost ./ addNoRecourseCost 100 * multExtCost ./ multNoRecourseCost], outliers = false, title = "(b) High capacity")

    outPlot = plot(plot1, plot2, fontfamily = "Times", xticks = (1:2, ["Additive" "Multiplicative"]), ylabel = "Relative cost [in %]", ylims = (70, 130), legend = false, fill = [:lightgrey :lightgrey], size = (800, 350), left_margin = 4mm)
    savefig(outPlot,"output/boxplotValueRecourse.pdf")

    return outPlot
end
