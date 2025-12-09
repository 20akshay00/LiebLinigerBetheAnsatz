function get_magnon_spectrum(γ; N=100)
    c = 1.0

    _, _, n_total, Q = get_ground_state(γ; N=N)
    k_fermi = π * n_total

    Ksym(k, q) = c / π * (1 / (c^2 + (k - q)^2) + 1 / (c^2 + (k + q)^2))
    solver = QuadratureSolver(gausslobatto(N))

    rho, xs, ws = solve(solver, Ksym, x -> 1 / (2π), 0.0, Q)
    rho_vals = rho.(xs)

    # solve for dressed energy epsilon
    # (I-K)ε₀ = k^2  and  (I-K)ε₁ = 1
    eps0_func, _, _ = solve(solver, Ksym, k -> k^2, 0.0, Q)
    eps1_func, _, _ = solve(solver, Ksym, k -> 1.0, 0.0, Q)

    # determine chemical potential μ such that ε(Q) = 0
    # ε(k) = ε₀(k) - μ * ε₁(k)
    μ = eps0_func(Q) / eps1_func(Q)
    eps_vals = eps0_func.(xs) .- μ .* eps1_func.(xs)

    # compute magnon Spectrum
    # scan spin rapidity Λ from 0 to ∞
    # we use a denser grid near 0 where curvature is high
    Λs = [range(0, 5.0, length=80); range(5.1, 50.0, length=20)]

    P_mag = Float64[]
    E_mag = Float64[]

    # kernel for magnon coupling (Lorentzian of width c)
    # K_coupling(x) = c/π * 1/(c^2 + x^2) 
    K_coupling(x) = (c / π) * (1 / (c^2 + x^2))

    theta(x) = 2 * atan(x / c)

    for Λ in Λs
        # energy: E = - ∫_{-Q}^Q K_coupling(k - Λ) ε(k) dk
        term_E = K_coupling.(xs .- Λ) .+ K_coupling.(xs .+ Λ)
        E = -dot(ws, term_E .* eps_vals)

        # momentum: P = k_F - ∫_{-Q}^Q θ(k - Λ) ρ(k) dk
        term_P = theta.(xs .- Λ) .- theta.(xs .+ Λ)
        P = k_fermi - dot(ws, term_P .* rho_vals)

        push!(P_mag, P)
        push!(E_mag, E)
    end

    return P_mag, E_mag, k_fermi
end