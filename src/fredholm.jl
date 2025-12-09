# A Computational Approach to the Fredholm
# Integral Equation of the Second Kind
# https://www.iaeng.org/publication/WCE2008/WCE2008_pp933-937.pdf

abstract type FredholmSolver end
struct QuadratureSolver <: FredholmSolver
    weights
    roots
    N::Int

    function QuadratureSolver((r, w),)
        @assert length(w) == length(r) "Weights and roots must have same length!"
        new(w, r, length(w))
    end
end

# f(x) - λ∫dy kernel(x, y) f(y) = g(x)
function solve(solver::QuadratureSolver, kernel, g, yi, yf)
    (; weights, roots, N) = solver

    # map roots and weights to the physical interval [yi, yf]
    xs, ws = rescale(roots, weights, yi, yf)

    K = kernel.(xs, transpose(xs))
    gs = g.(xs)
    fs = (I - (K .* transpose(ws))) \ gs

    return (x -> dot(ws .* fs, kernel.(x, xs)) + g(x)), xs, ws
end

### performs better for singularities or finite peaks/bumps around |x-y| → 0
struct ModifiedQuadratureSolver <: FredholmSolver
    weights
    roots
    N::Int

    function ModifiedQuadratureSolver((r, w),)
        @assert length(w) == length(r) "Weights and roots must have same length!"
        new(w, r, length(w))
    end
end

# f(x) - λ∫dy kernel(x, y) f(y) = g(x)
function solve(solver::ModifiedQuadratureSolver, kernel, g, kernel_traced, yi, yf)
    (; weights, roots, N) = solver

    # map roots and weights to the physical interval [yi, yf]
    xs, ws = rescale(roots, weights, yi, yf)

    K = kernel.(xs, transpose(xs))
    gs = g.(xs)

    Δ = Diagonal(K * ws .- kernel_traced.(xs))
    fs = (I + (Δ - K .* transpose(ws))) \ gs
    return (x -> (k = kernel.(x, xs); (g(x) + dot(ws .* fs, k)) / (1 + dot(ws, k) - kernel_traced(x)))), xs, ws
end
