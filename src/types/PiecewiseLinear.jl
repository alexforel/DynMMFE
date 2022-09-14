# This function was taken from the website "http://www.pkofod.com/", article "http://www.pkofod.com/2017/05/30/solving-a-simple-discrete-choice-model-using-gaussian-quadrature/"
# I adapted the function to extrapolate, that is evaluating outside of the bounds it was created with.
struct PiecewiseLinear
    nodes
    values
    slopes

    function PiecewiseLinear(nodes, f)
        slopes = Float64[]
        fn = f.(nodes)
        for i = 1:length(nodes)-1
            # node i and node i+1
            ni, nip1 = nodes[i], nodes[i+1]
            # f evaluated at the nodes
            fi, fip1 = fn[i], fn[i+1]
            # store slopes in each interval, so we don't have to recalculate them every time
            push!(slopes, (fip1-fi)/(nip1-ni))
        end
        # Construct the type
        new(nodes, fn, slopes)
    end


    function (p::PiecewiseLinear)(x)
        index_low = searchsortedlast(p.nodes, x)
        n = length(p.nodes)

        # Interpolation: returns values as piecewise linear approximation between two calculated nodes and values.
        if 0 < index_low < n
            # x is in in the nodes of p
            return p.values[index_low + 1] + (x - p.nodes[index_low + 1]) * p.slopes[index_low]
        end

        # Extrapolation: allows to evaluate the pwl at any point x between two nodes using the slopes
        # Note that, in the original version, the function is "truncated" in the sense that values of x outside of the breakpoints will return the nearest breakpoint value instead of extrapolating the out of scope value with the slope of the last linear segment.
        if index_low == n
            # x is bigger than the biggest node of p
            lastY = p.values[end]
            lastX = p.nodes[end]
            lastSlope = p.slopes[end]
            return lastY + lastSlope * (x - lastX)
        elseif index_low == 0
            # x is smaller than the smallest node of p
            lastY = p.values[1]
            lastX = p.nodes[1]
            lastSlope = p.slopes[1]
            return lastY + lastSlope * (x - lastX)
        end
    end
end
