using QuadGK

"""
    NonUniformLLProblem(c, V; N=nothing, μ=nothing, domain=(-10.0, 10.0))

Problem definition for a Lieb-Liniger gas in an external potential V(x) using the 
Local Density Approximation (LDA).

# Arguments
- `c::Float64`: Interaction strength.
- `V::Function`: External potential function V(x).
- `N::Union{Nothing,Float64}`: Total particle number for canonical ensemble.
- `μ::Union{Nothing,Float64}`: Global chemical potential for grand canonical ensemble.
- `domain::Tuple{Float64,Float64}`: Spatial region [x_min, x_max] for integration.
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
    μ::Float64
    rho_x::Function
    e_x::Function
    N_total::Float64
    E_total::Float64
end

"""
    solve(p::NonUniformLLProblem; eos_samples=40, maxiter=50, atol=1e-8)

Solves the non-uniform Lieb-Liniger problem by sampling the Equation of State (EoS) 
and integrating the local density over the spatial domain.

# Arguments
- `eos_samples`: Number of points used to sample the μ -> density/energy mapping.
- `maxiter`: Maximum iterations for the chemical potential bisection.
- `atol`: Absolute tolerance for bisection convergence.
"""
function solve(p::NonUniformLLProblem; eos_samples=40, maxiter=50, atol=1e-8)
    c, V = p.c, p.V
    x_min, x_max = p.domain

    # determine upper bound for chemical potential sampling
    μ_high = !isnothing(p.μ) ? p.μ : (isnothing(p.N) ? 10.0 : p.N * 0.5)

    # generate equation of state (EoS) lookup tables
    μ_grid = collect(range(1e-6, μ_high * 2.0, length=eos_samples))
    rho_eos = zeros(eos_samples)
    e_eos = zeros(eos_samples)

    # solve uniform system for each μ in the grid
    Threads.@threads for i in 1:eos_samples
        st = solve(InfiniteLLProblem(c=c, μ=μ_grid[i]))
        rho_eos[i] = average_particle_density(st)
        e_eos[i] = energy_density(st)
    end

    # local density and energy mapping via linear interpolation
    get_rho_from_mu(μ_loc) = μ_loc <= 0 ? 0.0 : _lerp(μ_loc, μ_grid, rho_eos)
    get_e_from_mu(μ_loc) = μ_loc <= 0 ? 0.0 : _lerp(μ_loc, μ_grid, e_eos)

    # determine global chemical potential
    μ_global = 0.0
    if !isnothing(p.μ)
        μ_global = p.μ
    else
        # bisection to find μ such that ∫ rho(μ - V(x)) dx = N
        target_n = p.N
        calc_n(m) = quadgk(x -> get_rho_from_mu(m - V(x)), x_min, x_max)[1]

        # bracket the root
        ml, mh = 0.0, μ_high
        iter = 0
        while calc_n(mh) < target_n
            mh *= 2.0
            iter += 1
            (iter > maxiter) && error("Failed to bracket μ for target N.")
        end

        # bisection loop
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

    # generate profile functions
    rho_x(x) = get_rho_from_mu(μ_global - V(x))
    e_int_x(x) = get_e_from_mu(μ_global - V(x))

    # calculate total particle number and energy
    # energy includes internal and potential energy densities
    n_final, _ = quadgk(rho_x, x_min, x_max)
    e_total, _ = quadgk(x -> e_int_x(x) + V(x) * rho_x(x), x_min, x_max)

    return NonUniformLLState(p, μ_global, rho_x, e_int_x, n_final, e_total)
end

energy(s::NonUniformLLState) = s.E_total
average_particle_density(s::NonUniformLLState) = s.N_total / (s.prob.domain[2] - s.prob.domain[1])
particle_density(s::NonUniformLLState) = s.rho_x
domain(s::NonUniformLLState) = s.prob.domain
quasimomentum_distribution(s::NonUniformLLState) = throw(DomainError("LDA does not support quasimomentum distributions."))
fermi_quasimomentum(s::NonUniformLLState) = throw(DomainError("LDA does not support fermi quasimomentum calculation."))