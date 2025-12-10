# solves f(x) - λ ∫([K(x,y) + K(x,-y)] f(y), 0, Q) dy = g(x)
# assumes that ρ(k) = ρ(-k) always
function solve_ll_distribution(c, Q; N=100, quadrature_rule=gausslobatto)
    # original kernel: 2c / (c^2 + (k-q)^2)
    # symmetrized:     K(k, q) + K(k, -q)
    kernel(k, q) = c / π * (1 / (c^2 + (k - q)^2) + 1 / (c^2 + (k + q)^2))
    kernel_traced(k) = 1 / π * (atan((k + Q) / c) - atan((k - Q) / c))

    rho, xs, ws = solve(
        ModifiedQuadratureSolver(quadrature_rule(N)),
        kernel,
        (x) -> 1. / (2π),
        kernel_traced,
        0.,
        Q
    )

    # density n = 2 * ∫(ρ(k), 0, Q) dk
    particle_density = 2 * dot(ws, rho.(xs))

    # energy E/L = 2 * ∫(k^2 ρ(k), 0, Q) dk
    energy_density = 2 * dot(ws, (xs .^ 2) .* rho.(xs))

    return rho, particle_density, energy_density
end

function get_ground_state(; γ=nothing, μ=nothing, c=1.0, Ql=1e-8, Qh=1., maxiter=100, kwargs...)
    if !isnothing(γ) && isnothing(μ)
        # particle density
        target = c / γ
        metric_func = Q -> solve_ll_distribution(c, Q; kwargs...)[2]
    elseif !isnothing(μ) && isnothing(γ)
        # chemical potential
        target = μ
        metric_func = Q -> compute_chemical_potential(c, Q; kwargs...)
    else
        error("Cannot specify both γ (canonical) and μ (grand canonical) at once!")
    end

    function residual(Q)
        (Q <= 1e-9) && return -target
        return metric_func(Q) - target
    end

    iter = 0
    while residual(Qh) < 0
        Qh *= 2.0
        iter += 1
        (iter > maxiter) && error("Could not bracket root (Target too high?)")
    end

    Q = find_zero(residual, (Ql, Qh), Bisection())

    rho, n, e = solve_ll_distribution(c, Q; kwargs...)
    return rho, e, n, Q
end

function compute_dressed_energy(c, Q; N=100, quadrature_rule=gausslobatto)
    kernel(k, q) = c / π * (1 / (c^2 + (k - q)^2) + 1 / (c^2 + (k + q)^2))
    kernel_traced(k) = 1 / π * (atan((k + Q) / c) - atan((k - Q) / c))

    solver = ModifiedQuadratureSolver(quadrature_rule(N))

    # ε(k) - ∫ K ε = k^2 - μ
    # solve auxiliary equations: (I - K)ε₀ = k^2  and  (I - K)ε₁ = 1
    eps0, _, _ = solve(solver, kernel, k -> k^2, kernel_traced, 0., Q)
    eps1, _, _ = solve(solver, kernel, k -> 1.0, kernel_traced, 0., Q)

    # enforce ε(Q) = 0 to find μ
    μ = eps0(Q) / eps1(Q)

    return (k) -> eps0(k) - μ * eps1(k), μ
end

function compute_chemical_potential(c, Q; N=100, quadrature_rule=gausslobatto, kwargs...)
    _, μ = compute_dressed_energy(c, Q; kwargs...)
    return μ
end

function get_particle_hole_spectrum(γ, c=1.; quadrature_rule=gausslobatto, N=100, num_points=100, kwargs...)
    rho_gs, _, _, Q = get_ground_state(γ=γ, c=c, kwargs...)

    ε, _ = compute_dressed_energy(c, Q, kwargs...)

    # dressed momentum P(k) = k + ∫ θ(k-q)ρ(q)dq
    xs, ws = rescale(quadrature_rule(N)..., 0., Q)
    θ(x) = 2 * atan(x / c)
    P(k) = k + dot(ws, (θ.(k .- xs) .+ θ.(k .+ xs)) .* rho_gs.(xs))

    kf = P(Q) # Fermi momentum

    # grid for half the Fermi sea [0, Q]
    k_h = range(0, Q, length=num_points)
    P_vals = P.(k_h)
    E_vals = -ε.(k_h) # energy is positive for excitations

    # #1 removing particle from right side (Q -> 0)
    # Momentum p = P(Q) - P(k)
    p1 = kf .- P_vals

    # #2 removing particle from left side (0 -> -Q)
    # Momentum p = P(Q) - P(-k) = P(Q) + P(k)
    p2 = kf .+ P_vals

    # full range 0 -> 2kF
    p_h = vcat(reverse(p1), p2)
    e_h = vcat(reverse(E_vals), E_vals)

    # Type II (Particle) branch (k > Q)
    k_p = range(Q, 3 * Q, length=N)
    p_p = P.(k_p) .- kf
    e_p = ε.(k_p)

    return p_h, e_h, p_p, e_p, kf
end