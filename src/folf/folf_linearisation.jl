function folf_linearisation(initInv, cumulatedDemandDistributionMatrix, nbSegments, capacity, nbProducts)
    ## folf_linearisation(initInv, cumulatedDemandDistributionMatrix, nbSegments, capacity, nbProducts)
    # Caclulate the slope of the first-order loss function and its complementary over each segment. Segments are deduced from the initial inventory and avaialble capacity.
    # Input:
    #   - initInv, cumulatedDemandDistributionMatrix, nbSegments, capacity, nbProducts
    # Output:
    #   - uRange: Matrix[nbProducts, T, nbSegments], value (of the cumulated production) at which the value of the folf and c-folf should be calculated
    #   - folfOrigin: Matrix[nbProducts, T], value of folf at first breakpoint
    #   - folfSlope: Matrix[nbProducts, T, nbSegments], slope of folf between two breakpoints
    #   - complFolfOrigin Matrix[nbProducts, T], value of c-folf at first breakpoint
    #   - complFolfSlope: Matrix[nbProducts, T, nbSegments], slope of c-folf between two breakpoints

    if VERBOSE == 1
        println("- Starting linearisation of FOLF -")
    end
    ## Determine segments over which to linearise
    uRange = zeros(nbProducts, T, nbSegments + 1)
    # Segments are determines over the whole feasible domain of inventory positions: from starting inventory to the max position with full production
    for k in 1:nbProducts
        for t in 1:T
            lowerLinBound = initInv[k,1]
            upperLinBound = initInv[k,1] + capacity * t
            uRange[k, t, :] = collect(range(lowerLinBound, stop = upperLinBound,
                                            length = (nbSegments+1)))
        end
    end

    # Calculate values of FOLF at selected points
    folfOrigin, folfValueMatrix, complFolfOrigin, complFolfValueMatrix = evaluate_folf_and_complementary(uRange, cumulatedDemandDistributionMatrix, nbSegments, capacity, nbProducts)

    # Linearise the first order loss function of all products, all time periods
    folfSlope, complFolfSlope = calculate_slope_of_folf_and_complementary(uRange, folfValueMatrix, complFolfValueMatrix, nbSegments, nbProducts)

    return uRange, folfOrigin, folfSlope, complFolfOrigin, complFolfSlope
end

function evaluate_folf_and_complementary(uRange, cumulatedDemandDistributionMatrix, nbSegments, capacity, nbProducts)
    folfValueMatrix = zeros(nbProducts, T, nbSegments + 1)
    complFolfValueMatrix = zeros(nbProducts, T, nbSegments + 1)
    # Calculate values of FOLF at selected points
    for k = 1:nbProducts
        for t = 1:T
            for l in 1:(nbSegments+1)
                folfValueMatrix[k, t, l] = first_order_loss_function(uRange[k, t, l], cumulatedDemandDistributionMatrix[k, t], capacity)
                complFolfValueMatrix[k, t, l] = folfValueMatrix[k, t, l] + uRange[k, t, l] - mean(cumulatedDemandDistributionMatrix[k,t])
            end
        end
    end

    # Origin = value at first breakpoint
    folfOrigin = folfValueMatrix[:, :, 1]
    complFolfOrigin = complFolfValueMatrix[:, :, 1]

    return folfOrigin, folfValueMatrix, complFolfOrigin, complFolfValueMatrix
end

function calculate_slope_of_folf_and_complementary(uRange, folfValueMatrix, complFolfValueMatrix, nbSegments, nbProducts)
    folfSlope = zeros(nbProducts, T, nbSegments)
    complFolfSlope = zeros(nbProducts, T, nbSegments)
    for k = 1:nbProducts
        for t = 1:T
            # Check that uRange is not only zeros
            if !(prod(uRange[k, t, :] .== 0))
                # Slope of FOLF
                folfSlope[k, t, :] = (folfValueMatrix[k,t, 2:end] - folfValueMatrix[k,t, 1:(end-1)]) ./ (uRange[k,t, 2:end] - uRange[k,t, 1:(end-1)])
                # Slope of c-FOLF
                complFolfSlope[k, t, :] = (complFolfValueMatrix[k,t, 2:end] - complFolfValueMatrix[k,t, 1:(end-1)]) ./ (uRange[k,t, 2:end] - uRange[k,t, 1:(end-1)])
            end
        end
    end

    return folfSlope, complFolfSlope
end



function free_returns_folf_linearisation(inventoryState, cumulatedDemandDistributionMatrix, nbSegments, capacity, nbProds)
    ## Determine segments over which to linearise
    vRange = zeros(nbProds, nbSegments + 1)
    # Segments are determines over the whole feasible domain of inventory positions: from starting inventory to the max position with full production
    for k in 1:nbProds
        # Note that the breakpoints are ordered in decreasing fashion
        lowerLinBound = max(inventoryState[k,1], 0)
        upperLinBound = 0
        vRange[k, :] = collect(range(lowerLinBound, stop = upperLinBound, length = (nbSegments+1)))
    end

    # Calculate values of FOLF at selected points
    folfValueMatrix = zeros(nbProds, nbSegments + 1)
    complFolfValueMatrix = zeros(nbProds, nbSegments + 1)
    # Calculate values of FOLF at selected points
    for k = 1:nbProds
        if !(prod(vRange[k, :] .== 0))
            for l in 1:(nbSegments+1)
                folfValueMatrix[k, l] = first_order_loss_function(vRange[k, l], cumulatedDemandDistributionMatrix[k, 1], capacity)
                complFolfValueMatrix[k, l] = folfValueMatrix[k, l] + vRange[k, l] - mean(cumulatedDemandDistributionMatrix[k,1])
            end
        end
    end

    # Origin = value at first breakpoint
    folfOrigin = folfValueMatrix[:, 1]
    complFolfOrigin = complFolfValueMatrix[:, 1]

    # Linearise the first order loss function of all products, all time periods
    folfSlope = zeros(nbProds, nbSegments)
    complFolfSlope = zeros(nbProds, nbSegments)
    for k = 1:nbProds
        if !(prod(vRange[k, :] .== 0))
            # Slope of FOLF
            folfSlope[k, :] = (folfValueMatrix[k, 2:end] - folfValueMatrix[k, 1:(end-1)]) ./ (vRange[k, 2:end] - vRange[k, 1:(end-1)])
            # Slope of c-FOLF
            complFolfSlope[k, :] = (complFolfValueMatrix[k, 2:end] - complFolfValueMatrix[k, 1:(end-1)]) ./ (vRange[k, 2:end] - vRange[k, 1:(end-1)])
        end
    end

    return vRange, folfOrigin, folfSlope, complFolfOrigin, complFolfSlope
end

function check_return_linearization_consistency(vRange, returnFolfOrigin, folfOrigin, returnComplFolfOrigin, complFolfOrigin, nbProds)
	# Quick consistency check
	for k in 1:nbProds
		if !(prod(vRange[k, :] .== 0))
			if !(returnFolfOrigin[k] == folfOrigin[k,1])
				print("Warning: FOLF origin different for production and returns for product k=", k, " .")
			end
			if !(returnComplFolfOrigin[k] == complFolfOrigin[k,1])
				print("Warning: FOLF origin different for production and returns for product k=", k, " .")
			end
		end
	end
end
