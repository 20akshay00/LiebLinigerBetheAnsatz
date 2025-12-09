using DrWatson
@quickactivate "LiebLinigerBetheAnsatz"
using LiebLinigerBetheAnsatz
using LinearAlgebra, Plots, LaTeXStrings

begin
    # TG limit
    rho, e, n, Q = get_ground_state(1e8, N=100)
    println(norm(e / n^3 - π^2 / 3), " ", Q, " ", rho(0))
    plot(rho, range(-Q, Q, 100), ylim=[0, 4.], lab="", lw=4)
end

let
    γ = 10
    N = 100
    quadrature_rule = midpoint_quadrature
    rho, e, n, Q = get_ground_state(γ, N=N, quadrature_rule=quadrature_rule)
    p_h, E_h, p_p, E_p = get_excitation_spectrum(γ, N=N, quadrature_rule=quadrature_rule)

    # what is this factor of 4??
    plot(p_p ./ n, E_p ./ n^2, label="Type II (Particles)", lw=2, xtick=pitick(0, 2π, 1, mode=:latex))
    Eph(p) = sqrt(4 * γ * n^2 * p^2 + p^4)
    plot!(p_p ./ n, Eph.(p_p) ./ n^2, ls=:dash, c=:black, lw=2, lab="BdG")

    c = sqrt(4 * γ * n^2)
    Esol(v) = (4 / 3) * n * c * (1 - (v / c)^2)^(3 / 2)
    psol(v) = 2 * n * (acos(v / c) - (v / c) * sqrt(1 - (v / c)^2))
    vs = range(-c, c, 100)

    plot!(p_h ./ n, E_h ./ n^2, label="Type I (Holes)", lw=2, title="Lieb-Liniger Spectrum (γ=10)")
    plot!(psol.(vs) ./ n, Esol.(vs) ./ n^2, c=:black, ls=:dash, lw=2, lab="Soliton")

    xlabel!("Momentum " * L"p/ρ")
    ylabel!("Energy " * L"\epsilon/ρ^2")
end

