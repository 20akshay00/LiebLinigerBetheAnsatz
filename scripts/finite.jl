using LiebLinigerBetheAnsatz
using LinearAlgebra, Plots, LaTeXStrings

begin # single parameter
    @time state = solve(FiniteLLProblem(1., 1e-7, N=2, bc=:hardwall))
    rho_k = quasimomentum_distribution(state)
    e = energy_density(state)
    n = average_particle_density(state)
    Q = fermi_quasimomentum(state)
    println(energy(state))
end

begin # convergence to thermodynamic limit
    ρ, c = 10., 10.
    bc = :periodic
    Ns = 100:100:600
    energy_densities = zeros(length(Ns))

    for idx in eachindex(Ns)
        N = Ns[idx]
        @time state = solve(FiniteLLProblem(N / ρ, c, N=N, bc=bc))
        energy_densities[idx] = energy_density(state)
    end

    @time state_tdl = solve(InfiniteLLProblem(γ=c / ρ, c=c))
end;

begin
    scatter(Ns, energy_densities, lab="bc=$(string(bc))")
    hline!([energy_density(state_tdl)], ls=:dash, c=:black, lw=2, lab="")
    plot!(title="Convergence to TDL | ρ = $ρ | c = $c", xlabel="Particle number", ylabel="Energy density", framestyle=:box)
end

begin # log-log
    scatter(Ns, abs.(energy_densities .- energy_densities[end]) .+ 1e-8, lab="bc=$(string(bc))", yscale=:log10)
    plot!(title="Convergence to TDL | ρ = $ρ | c = $c", xlabel="Particle number", ylabel="Energy density", framestyle=:box)
end



