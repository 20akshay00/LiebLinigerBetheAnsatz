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

function get_ground_state(γ; c=1.0, kwargs...)
    particle_density_target = c / γ

    function residual(Q)
        (Q <= 1e-9) && return -particle_density_target
        _, particle_density_curr, _ = solve_ll_distribution(c, Q; kwargs...)
        return particle_density_curr - particle_density_target
    end

    Q_low = 1e-8 # n ~ 0 => residual(Q) < 0
    Q_high = 1.0 # find such that residual(Q) > 0

    # expand upper bound until density is high enough
    iter = 0
    while residual(Q_high) < 0
        Q_high *= 2.0
        iter += 1
        if iter > 50
            error("Could not bracket root (Target density too high?)")
        end
    end

    Q = find_zero(residual, (Q_low, Q_high), Bisection())

    rho, particle_density_final, energy_density_final = solve_ll_distribution(c, Q; kwargs...)

    return rho, energy_density_final, particle_density_final, Q
end

function get_ground_state_gce(μ, c; N=100, kwargs...)
    function residual(Q)
        (Q <= 1e-9) && return -μ
        return compute_chemical_potential(c, Q; N=N, kwargs...) - μ
    end

    Q_low, Q_high = 1e-6, 10.0
    while residual(Q_high) < 0
        Q_high *= 2.0
    end

    Q = find_zero(residual, (Q_low, Q_high), Bisection())

    rho, particle_density_final, energy_density_final = solve_ll_distribution(c, Q; N=N, kwargs...)
    return rho, energy_density_final, particle_density_final, Q
end

function compute_dressed_energy(c, Q; N=100, quadrature_rule=gausslobatto)
    K(k, q) = c / π * (1 / (c^2 + (k - q)^2) + 1 / (c^2 + (k + q)^2))

    solver = QuadratureSolver(quadrature_rule(N))

    # ε(k) - ∫ K ε = k^2 - μ
    # solve auxiliary equations: (I - K)ε₀ = k^2  and  (I - K)ε₁ = 1
    eps0, _, _ = solve(solver, K, k -> k^2, 0., Q)
    eps1, _, _ = solve(solver, K, k -> 1.0, 0., Q)

    # enforce ε(Q) = 0 to find μ
    μ = eps0(Q) / eps1(Q)

    return (k) -> eps0(k) - μ * eps1(k), μ
end

function compute_chemical_potential(c, Q; N=100, quadrature_rule=gausslobatto, kwargs...)
    _, μ = compute_dressed_energy(c, Q; kwargs...)
    return μ
end

function get_excitation_spectrum(γ, c=1.; quadrature_rule=gausslobatto, N=100, kwargs...)
    rho_gs, _, _, Q = get_ground_state(γ, c=c; kwargs...)

    ε, _ = compute_dressed_energy(c, Q, kwargs...)

    # dressed momentum P(k) = k + ∫ θ(k-q)ρ(q)dq
    xs, ws = rescale(quadrature_rule(N)..., 0., Q)
    θ(x) = 2 * atan(x / c)
    P(k) = k + dot(ws, (θ.(k .- xs) .+ θ.(k .+ xs)) .* rho_gs.(xs))

    # dispersion Curves
    k_fermi = P(Q)

    # grid for half the Fermi sea [0, Q]
    k_h = range(0, Q, length=N)
    P_vals = P.(k_h)
    E_vals = -ε.(k_h) # energy is positive for excitations

    # #1 removing particle from right side (Q -> 0)
    # Momentum p = P(Q) - P(k)
    p1 = k_fermi .- P_vals

    # #2 removing particle from left side (0 -> -Q)
    # Momentum p = P(Q) - P(-k) = P(Q) + P(k)
    p2 = k_fermi .+ P_vals

    # full range 0 -> 2kF
    p_h = vcat(reverse(p1), p2)
    e_h = vcat(reverse(E_vals), E_vals)

    # Type II (Particle) branch (k > Q)
    k_p = range(Q, 3 * Q, length=N)
    p_p = P.(k_p) .- k_fermi
    e_p = ε.(k_p)

    return p_h, e_h, p_p, e_p, k_fermi
end