## common quadrature rules
# standardized to [-1, 1]
function reimann_quadrature(N)
    xs = range(-1, 1, N)
    return xs, fill(step(xs), N)
end

function midpoint_quadrature(N)
    dx = 2 / N
    xs = range(-1 + dx / 2, 1 - dx / 2, N)
    return xs, fill(dx, N)
end

function trapezoidal_quadrature(N)
    xs = range(-1, 1, N)
    dx = step(xs)
    ws = fill(dx, N)
    ws[1] /= 2
    ws[end] /= 2
    return xs, ws
end

function simpson_quadrature(N)
    xs = range(-1, 1, N)
    dx = step(xs)
    ws = zeros(N)

    i = 1
    while i <= N - 1
        if N - i + 1 >= 4 && (N - i + 1) % 2 == 0  # use 3/8 rule for 3 intervals if needed
            ws[i] += 3dx / 8
            ws[i+1] += 9dx / 8
            ws[i+2] += 9dx / 8
            ws[i+3] += 3dx / 8
            i += 3
        else  # use 1/3 rule for 2 intervals
            ws[i] += dx / 3
            ws[i+1] += 4dx / 3
            ws[i+2] += dx / 3
            i += 2
        end
    end

    return xs, ws
end

## axis labels as multiples of pi
function piticklabel(x::Rational, ::Val{:latex})
    iszero(x) && return L"0"
    S = x < 0 ? "-" : ""
    n, d = abs(numerator(x)), denominator(x)
    N = n == 1 ? "" : repr(n)
    d == 1 && return L"%$S%$N\pi"
    L"%$S\frac{%$N\pi}{%$d}"
end

function pitick(start, stop, denom; mode=:text)
    # Compute lower index
    a = Int(floor(start / (π / denom)))
    # Compute upper index
    b = Int(ceil(stop / (π / denom)))
    tick = range(a * π / denom, b * π / denom; step=π / denom)
    ticklabel = piticklabel.((a:b) .// denom, Val(mode))
    tick, ticklabel
end