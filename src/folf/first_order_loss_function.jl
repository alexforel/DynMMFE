function first_order_loss_function(x, probaDistribution, capacity)
    ## first_order_loss_function(x, probaDistribution)
    # Evaluate the first-order loss function of an input probaDistribution at value x.
    # Input:
    #   - x, probaDistribution
    #   - capacity: scalar, used to calculate upper bound of integral
    # Output:
    #   - folfValue

    REL_TOL_INTEGRAL = 1e-9

    if typeof(probaDistribution) == Normal{Float64}
        # FOLF for normal distribution can be calculated using FOLF standard normal
        tempMean = mean(probaDistribution)
        tempStd  = std(probaDistribution)
        folfValue = tempStd * pwlaStandardFolf((x - tempMean) / tempStd)
    elseif typeof(probaDistribution) == LogNormal{Float64}
        if x < 0
            residual = abs(x)
            xMin = 0.0
        else
            residual = 0.0
            xMin = x
        end
        integrand = simpson_integration(t -> 1 - cdf(probaDistribution, t), xMin, 5 * capacity, 2000)
        folfValue = integrand + residual
    elseif typeof(probaDistribution) == Truncated{Normal{Float64},Continuous,Float64}
        integrand = simpson_integration(t -> 1 - cdf(probaDistribution, t), x, 5 * capacity, 2000)
        folfValue = integrand
    end

    return folfValue
end

function first_order_loss_function(x, probaDistribution)
    first_order_loss_function(x, probaDistribution, 2 * mean(probaDistribution))
end
