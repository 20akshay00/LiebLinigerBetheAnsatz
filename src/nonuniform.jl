using QuadGK
using LinearAlgebra

"""
    NonUniformLLProblem(; c, V, N=nothing, μ=nothing, domain=(-1.0, 1.0))

Problem definition for a Lieb-Liniger gas in an external potential V(x) using the 
Local Density Approximation (LDA) with a spectral approximation of the equation of state (EoS) ρ(μ) and e(μ).

# Arguments
- `c`: Interaction strength.
- `V`: External potential function V(x).
- `N`: Total particle number.
- `μ`: Global chemical potential.
- `domain`: The spatial region [x_min, x_max].
"""
struct NonUniformLLProblem <: LLProblem
    c::Float64
    V::Function
    N::Union{Nothing,Float64}
    μ::Union{Nothing,Float64}
    domain::Tuple{Float64,Float64}

    function NonUniformLLProblem(; c, V, N=nothing, μ=nothing, domain=(-1.0, 1.0))
        if isnothing(N) == isnothing(μ)
            throw(ArgumentError("Must specify exactly one of N or μ."))
        end
        new(c, V, N, μ, domain)
    end
end

struct NonUniformLLState <: LLState
    prob::NonUniformLLProblem
    μ_global::Float64
    rho_x::Function
    e_x::Function
    N_total::Float64
    E_total::Float64
end

## Helper for Chebyshev functions

# generates M chebyshev nodes in the interval [-1, 1]
_cheb_nodes(M) = [cos(pi * (2k - 1) / (2M)) for k in 1:M]

# calculates chebyshev coefficients using a discrete cosine transform
function _cheb_coeffs(ys)
    M = length(ys)
    coeffs = zeros(M)
    for j in 0:M-1
        # normalization factor for chebyshev polynomials
        norm = (j == 0) ? 1 / M : 2 / M
        for k in 1:M
            coeffs[j+1] += ys[k] * cos(pi * j * (2k - 1) / (2M))
        end
        coeffs[j+1] *= norm
    end
    return coeffs
end

# evaluates the chebyshev expansion using clenshaw's recurrence algorithm
function _cheb_eval(x, coeffs, x_min, x_max)
    # map x from [x_min, x_max] to t in [-1, 1]

    t = (2 * x - (x_max + x_min)) / (x_max - x_min)
    t = clamp(t, -1.0, 1.0)

    M = length(coeffs)
    d1 = 0.0
    d2 = 0.0
    for j in M:-1:2
        d1, d2 = 2 * t * d1 - d2 + coeffs[j], d1
    end
    return t * d1 - d2 + coeffs[1]
end

"""
    solve(p::NonUniformLLProblem; n_spectral=20, maxiter=50, atol=1e-10)

Solves the NonUniform problem by constructing a spectral approximation of the EoS

# Arguments
- `n_spectral`: Number of Chebyshev nodes.
- `maxiter`: Maximum iterations for the μ_global bisection.
- `atol`: Absolute tolerance for μ_global.
"""
function solve(p::NonUniformLLProblem; n_spectral=25, maxiter=50, atol=1e-10)
    c, V = p.c, p.V
    x_min, x_max = p.domain

    # determine the range for our EoS surrogate [0, μ_high]
    μ_high = !isnothing(p.μ) ? p.μ : (isnothing(p.N) ? 10.0 : p.N * 0.5)
    μ_range = (1e-7, μ_high * 2.5)

    ## (1) Solve TDL at spectral nodes
    # map chebyshev nodes from [-1, 1] to the μ_range
    t_nodes = _cheb_nodes(n_spectral)
    μ_nodes = [0.5 * (μ_range[2] + μ_range[1]) + 0.5 * (μ_range[2] - μ_range[1]) * t for t in t_nodes]

    rho_samples = zeros(n_spectral)
    e_samples = zeros(n_spectral)

    Threads.@threads for i in 1:n_spectral
        st = solve(InfiniteLLProblem(c=c, μ=μ_nodes[i]))
        rho_samples[i] = average_particle_density(st)
        e_samples[i] = energy_density(st)
    end

    # compute spectral coefficients for the EoS
    rho_coeffs = _cheb_coeffs(rho_samples)
    e_coeffs = _cheb_coeffs(e_samples)

    # spectral approximations
    # we return 0.0 if μ <= 0 (the vacuum region)
    get_rho(μ_loc) = μ_loc <= μ_range[1] ? 0.0 : _cheb_eval(μ_loc, rho_coeffs, μ_range[1], μ_range[2])
    get_e(μ_loc) = μ_loc <= μ_range[1] ? 0.0 : _cheb_eval(μ_loc, e_coeffs, μ_range[1], μ_range[2])

    ## (2) Solve for μ_global
    μ_global = 0.0
    if !isnothing(p.μ)
        μ_global = p.μ
    else
        target_n = p.N
        # calculate total N using the spectral surrogate
        calc_n(m) = quadgk(x -> get_rho(m - V(x)), x_min, x_max)[1]

        # bracket and bisect
        ml, mh = 0.0, μ_high
        iter = 0
        while calc_n(mh) < target_n
            mh *= 2.0
            iter += 1
            (iter > maxiter) && error("Failed to bracket μ_global.")
        end

        for _ in 1:maxiter
            mid = (ml + mh) / 2
            if calc_n(mid) < target_n
                ml = mid
            else
                mh = mid
            end
            abs(mh - ml) < atol && break
        end
        μ_global = (ml + mh) / 2
    end

    # (3) Construct profiles and integrate totals
    rho_x(x) = get_rho(μ_global - V(x))
    e_int_x(x) = get_e(μ_global - V(x))

    n_final, _ = quadgk(rho_x, x_min, x_max)
    e_total, _ = quadgk(x -> e_int_x(x) + V(x) * rho_x(x), x_min, x_max)

    return NonUniformLLState(p, μ_global, rho_x, e_int_x, n_final, e_total)
end

# api
energy(s::NonUniformLLState) = s.E_total
energy_density(s::NonUniformLLState) = s.E_total / (s.prob.domain[2] - s.prob.domain[1])
average_particle_density(s::NonUniformLLState) = s.N_total / (s.prob.domain[2] - s.prob.domain[1])
particle_density(s::NonUniformLLState) = s.rho_x
domain(s::NonUniformLLState) = s.prob.domain