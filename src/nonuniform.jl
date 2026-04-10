"""
    precompute_chebyshev_eos(c, μ_max; atol=1e-11, max_nodes=512)

Adaptively constructs Chebyshev interpolations for ρ(μ) and e(μ) (i.e, the Equation of State).

# Arguments
- `c`: Interaction strength.
- `μ_max`: Maximum chemical potential for the window.
- `atol`: Tolerance for the last Chebyshev coefficient.
- `max_nodes`: Limit on the number of nodes.
"""
function _construct_chebyshev_eos(c::Float64, μ_max::Float64; atol=1e-5, max_nodes=512)
    μ_min, n = 1e-8, 16
    while n <= max_nodes
        nodes = chebpoints(n, μ_min, μ_max)
        vals = [(st = solve(InfiniteLLProblem(c=c, μ=m));
        (average_particle_density(st), energy_density(st))) for m in nodes]

        rho_p = chebinterp(first.(vals), μ_min, μ_max)
        e_p = chebinterp(last.(vals), μ_min, μ_max)

        if abs(rho_p.coefs[end]) < atol && abs(e_p.coefs[end]) < atol
            rho_fn = μ -> μ <= μ_min ? 0.0 : rho_p(μ)
            e_fn = μ -> μ <= μ_min ? 0.0 : e_p(μ)

            @info "Chebyshev EoS converged to atol=$atol with $n nodes."
            return (rho=rho_fn, e=e_fn)
        end
        n *= 2
    end
    error("Chebyshev EoS failed to converge within $max_nodes nodes.")
end

"""
    NonUniformLLProblem(; c, V, N=nothing, μ=nothing, domain=(-1.0, 1.0), rho_eos=nothing, e_eos=nothing)

Problem definition for a Lieb-Liniger gas in an external potential V(x) using LDA.

# Arguments
- `c`: Interaction strength.
- `V`: External potential function V(x).
- `N`: Total particle number.
- `μ`: Global chemical potential.
- `domain`: The spatial region [x_min, x_max].
- `rho_eos`: A callable `f(μ)` for particle density.
- `e_eos`: A callable `f(μ)` for internal energy density.
"""
struct NonUniformLLProblem{FR,FE} <: LLProblem
    c::Float64
    V::Function
    N::Union{Nothing,Float64}
    μ::Union{Nothing,Float64}
    domain::Tuple{Float64,Float64}
    rho_eos::FR
    e_eos::FE

    function NonUniformLLProblem(; c, V, N=nothing, μ=nothing, domain=(-1.0, 1.0), rho_eos=nothing, e_eos=nothing)
        if isnothing(N) == isnothing(μ)
            throw(ArgumentError("Must specify exactly one of N or μ."))
        end

        if isnothing(rho_eos) || isnothing(e_eos)
            μ_est = !isnothing(μ) ? μ : (isnothing(N) ? 10.0 : N * 0.5)
            eos_pair = _construct_chebyshev_eos(c, μ_est * 2.5)
            rho_eos, e_eos = eos_pair.rho, eos_pair.e
        end

        new{typeof(rho_eos),typeof(e_eos)}(c, V, N, μ, domain, rho_eos, e_eos)
    end
end

struct NonUniformLLState{P<:NonUniformLLProblem} <: LLState
    prob::P
    μ0::Float64
    rho_x::Function
    e_x::Function
    N_total::Float64
    E_total::Float64
end

"""
    solve(p::NonUniformLLProblem; maxiter=50, atol=1e-10)

Solves the NonUniform problem via LDA integration.

# Arguments
- `p`: The problem definition.
- `maxiter`: Maximum iterations for bisection.
- `atol`: Absolute tolerance for convergence.
"""
function solve(p::NonUniformLLProblem; maxiter=50, atol=1e-10)
    x_min, x_max = p.domain

    μ0 = if !isnothing(p.μ)
        p.μ
    else
        calc_n(m) = quadgk(x -> p.rho_eos(m - p.V(x)), x_min, x_max)[1]

        ml, mh = 0.0, 1.0
        iter = 0
        while calc_n(mh) < p.N
            mh *= 2.0
            (iter += 1) > maxiter && error("Failed to bracket μ0.")
        end

        for _ in 1:maxiter
            mid = (ml + mh) / 2
            calc_n(mid) < p.N ? ml = mid : mh = mid
            abs(mh - ml) < atol && break
        end
        (ml + mh) / 2
    end

    rho_x(x) = p.rho_eos(μ0 - p.V(x))
    e_int_x(x) = p.e_eos(μ0 - p.V(x))

    n_final, _ = quadgk(rho_x, x_min, x_max)

    e_total, _ = quadgk(x_min, x_max) do x
        vx = p.V(x)
        μ_loc = μ0 - vx
        return p.e_eos(μ_loc) + vx * p.rho_eos(μ_loc)
    end

    return NonUniformLLState(p, μ0, rho_x, e_int_x, n_final, e_total)
end

energy(s::NonUniformLLState) = s.E_total
energy_density(s::NonUniformLLState) = s.E_total / (s.prob.domain[2] - s.prob.domain[1])
particle_density(s::NonUniformLLState) = s.rho_x
average_particle_density(s::NonUniformLLState) = s.N_total / (s.prob.domain[2] - s.prob.domain[1])
domain(s::NonUniformLLState) = s.prob.domain
chemical_potential(s::NonUniformLLState) = s.μ0