function simpson_integration(f::Function, a::Number, b::Number, n::Integer)
    # This function was taken from:
    # https://mmas.github.io/simpson-integration-julia
    h = (b - a) / n
    s = f(a + h / 2)
    for i in 1:(n-1)
        s += f(a + h * i + h / 2) + f(a + h * i) / 2
    end
    return h/6 * (f(a) + f(b) + 4*s)
end
