using LiebLinigerBetheAnsatz
using LinearAlgebra, Plots, LaTeXStrings

# helper to return ρ, ϵ, n, Q
get_observables(state) = (quasimomentum_distribution(state), energy_density(state), average_particle_density(state), fermi_quasimomentum(state))

begin
    # TG limit
    state = solve(InfiniteLLProblem(γ=1e7))
    rho, e, n, Q = get_observables(state)

    println(norm(e / n^3 - π^2 / 3), " ", Q, " ", rho(0))
    plot(rho, range(-Q, Q, 100), ylim=[0, 4.], lab="", lw=4)
end

let
    # excitations in dilute limit
    γ = 0.1
    N = 100
    quadrature_rule = midpoint_quadrature

    state = solve(InfiniteLLProblem(γ=γ), N=N, quadrature_rule=quadrature_rule)
    rho, e, n, Q = get_observables(state)
    p_h, E_h, p_p, E_p, = excitation_spectrum(state, N=N, quadrature_rule=quadrature_rule, num_points=20)

    plot()
    # what is this factor of 4??
    Eph(p) = sqrt(4 * γ * n^2 * p^2 + p^4)
    plot!(p_p ./ n, Eph.(p_p) ./ n^2, ls=:dash, c=:black, lw=1.5, lab="BdG")
    scatter!(p_p ./ n, E_p ./ n^2, label="Type II (Particles)", lw=2, xtick=pitick(0, 2π, 1, mode=:latex), c=1, ms=3)

    c = sqrt(4 * γ * n^2)
    Esol(v) = (4 / 3) * n * c * (1 - (v / c)^2)^(3 / 2)
    psol(v) = 2 * n * (acos(v / c) - (v / c) * sqrt(1 - (v / c)^2))
    vs = range(-c, c, 100)

    plot!(psol.(vs) ./ n, Esol.(vs) ./ n^2, c=:black, ls=:dash, lw=1.5, lab="Soliton")
    scatter!(p_h ./ n, E_h ./ n^2, label="Type I (Holes)", lw=2, title="Lieb-Liniger Spectrum (γ=$γ)", c=0, ms=3)

    xlabel!("Momentum " * L"p/ρ")
    ylabel!("Energy " * L"\epsilon/ρ^2")
    plot!(framestyle=:box)

    #savefig("./data/plots/particle-hole-gamma=$(γ)_c=1.png")
    # μ = compute_chemical_potential(c, Q)
    # println(μ, " ", 2 * n * c)
end

begin # speed of sound
    γs = 10 .^ (range(-2, 2, 20))
    vs = zeros(length(γs))
    for (idx, γ) in enumerate(γs)
        state = solve(InfiniteLLProblem(γ=γ))
        rho, e, n, Q = get_observables(state)
        vs[idx] = n / (2π * rho(Q)^2) / (2 * π * n)
    end
    plot(γs, vs, xscale=:log10, lw=2, lab="")
    plot!(ylims=[0, 1], ylabel=L"v_s/v_F", xlabel=L"\gamma")
end

begin # grand-canonical
    μ, c = 5., 10.
    state = solve(InfiniteLLProblem(μ=μ, c=c))
    _, e, n, _ = get_observables(state)

    println(e - μ * n)
    p_h, E_h, p_p, E_p = excitation_spectrum(state)

    plot(p_h ./ n, E_h ./ n^2, label="Type I (Holes)", lw=2, title="Lieb-Liniger Spectrum (γ=$γ)")
    plot!(p_p ./ n, E_p ./ n^2, label="Type II (Particles)", lw=2, xtick=pitick(0, 2π, 1, mode=:latex))
    xlabel!("Momentum " * L"p/ρ")
    ylabel!("Energy " * L"\epsilon/ρ^2")
end
