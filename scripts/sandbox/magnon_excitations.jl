using LiebLinigerBetheAnsatz
using LinearAlgebra, Plots, LaTeXStrings

# For the two-species bosonic model with SU(2) symmetric interactions 
# (i.e, Yang-Gaudin model), the ground state breaks the symmetry and 
# becomes polarized. This means that we may simply solve the single
# component Lieb-Liniger and all SU(2) rotations are valid ground states. 
# In this script, we utilize this ground state to compute the magnon 
# excitation branch by introducing particles of the other spin species.

get_observables(state) = (quasimomentum_distribution(state), energy_density(state), average_particle_density(state), fermi_quasimomentum(state))

begin # compare magnon w/ holon
    γ, c = 0.1, 1.
    state = solve(InfiniteLLProblem(γ=γ, c=c))
    rho, e, n, Q = get_observables(state)

    p_h, E_h, p_p, E_p, = excitation_spectrum(state)
    p_m, E_m, kf = magnon_spectrum(state)

    plot(p_h ./ n, E_h ./ n^2, label="Type I (Holes)", lw=2, title="Lieb-Liniger Spectrum (γ=$γ)")
    plot!(p_p ./ n, E_p ./ n^2, label="Type II (Particles)", lw=2, xtick=pitick(0, 2π, 1, mode=:latex))
    plot!(p_m ./ n, E_m ./ n^2, label="Magnon", lw=2, xtick=pitick(0, 2π, 1, mode=:latex))

    println(E_h[50] / E_m[50])
    xlabel!("Momentum " * L"p/ρ")
    ylabel!("Energy " * L"\epsilon/ρ^2")
    # ylims!(0, 20)
end

begin # magnon spectrum
    γ = 1.
    c = 1.

    state = solve(InfiniteLLProblem(γ=γ, c=c))
    rho, e, n, Q = get_observables(state)
    p, e, kf = magnon_spectrum(state)

    plot(p ./ n, e ./ n^2, ylim=[-0.1, Inf], lw=2, lab="", xtick=pitick(0, 2π, 1, mode=:latex), xlabel="", title="Lieb-Liniger Spectrum (γ=$γ)")
    plot!(framestyle="box")
    xlabel!("Momentum " * L"p/ρ")
    ylabel!("Energy " * L"\epsilon/ρ^2")
    #savefig("./data/plots/magnon-gamma=$(γ)_c=$c.png")
end


let # magnon spectrum
    γs = 10 .^ range(-3, 4, 40)
    c = 1

    p = plot()
    res = []
    for (i, γ) in enumerate(γs)
        state = solve(InfiniteLLProblem(γ=γ, c=c))
        rho, e, n, Q = get_observables(state)
        p, e, kf = magnon_spectrum(state)

        push!(res, e[1] ./ n^2)
        #plot!(p ./ n, e ./ n^2, ylim=[-0.1, Inf], lw=2, xtick=pitick(0, 2π, 1, mode=:latex), xlabel="", lab="γ=$γ", color = palette(:roma, length(γs))[i])
    end

    scatter(γs, res, xscale=:log10, yscale=:log10)
    plot!(framestyle="box")
    xlabel!("LL parameter " * L"\gamma")
    ylabel!("Maximum magnon energy " * L"\epsilon/ρ^2")
    #savefig("./data/plots/magnon-gamma=$(γ)_c=$c.png")
end
