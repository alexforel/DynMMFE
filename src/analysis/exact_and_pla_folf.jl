function exact_and_pla_folf(inputDistribution:: Distribution)
    ## exact_and_pla_folf(inputDistribution:: Distribution)
    # Compare the exact first-order loss function and its piecewise-linear approximation. Plot functions and return values.
    # Input:
    #   - inputDistribution: Distribution, distribution of (cumulative) demand
    # Output:
    #   - xVect, trueFolf, trueComplfFolf, uRange, linFolf, linComplFolf

    # Domain bounds on inventory for which to calculate folf
    xMin = 50
    xMax = 150
    xVect = collect(range(xMin, stop = xMax, length =100))

    ## --- First-order loss function ---
    trueFolf = first_order_loss_function.(xVect, inputDistribution)
    trueComplfFolf = trueFolf + (xVect .- mean(inputDistribution))
    # Piecewise-Linearisation of FOLF
    uRange, folfOrigin, folfSlope, complFolfOrigin, complFolfSlope = folf_linearisation(xMin, [inputDistribution], nbSegments, xMax - xMin, 1) # provide breakproints and slopes of folf and compl-folf
    nbPoints = nbSegments + 1
    # Get linearised functions
    linFolf = zeros(nbPoints) ; linComplFolf = zeros(nbPoints)
    linFolf[1] = folfOrigin[1,1] ; linComplFolf[1] = complFolfOrigin[1,1]
    for p in 1:nbSegments
        linFolf[p+1] = linFolf[p] + (uRange[p+1]-uRange[p]) * folfSlope[p]
        linComplFolf[p+1] = linComplFolf[p] + (uRange[p+1]-uRange[p]) * complFolfSlope[p]
    end
    # Plot results
    p1 = plot(xVect, trueFolf)
    p1 = plot!(uRange[1,1,:], linFolf)
    p1 = plot!(xVect, trueComplfFolf)
    p1 = plot!(uRange[1,1,:], linComplFolf)
    display(p1)

    return xVect, trueFolf, trueComplfFolf, uRange, linFolf, linComplFolf
end
