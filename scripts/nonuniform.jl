using LiebLinigerBetheAnsatz
using LinearAlgebra, Plots, LaTeXStrings

begin
    V = x -> 0.5x^2
    L = 9.
    μ = 10.
    c = 20.
end

begin
    @time state = solve(NonUniformLLProblem(c=c, V=V, μ=μ, domain=(-L, L)), n_spectral=100)

    dens = particle_density(state)
    xs = range(domain(state)..., 100)
    plot(dens, xs, lab="E = $(round(energy(state), digits=3))", lw=2)
    plot!(xlabel="Position", ylabel="Single particle density", framestyle=:box)
end

begin
    e = energy(state) * 0.5
    n = sum(dens.(xs) * step(xs))
    f = e - μ * n

    println("-----------------")
    println("μ = $μ | c = $c")
    println("-----------------")
    println("Energy: $(round(e, digits=3))")
    println("Particle number: $(round(n, digits=3))")
    println("Free energy: $(round(f, digits=3))")
end

begin
    # convergence with chebyshev nodes
    ns = 25:5:125
    res = zeros(length(ns))
    for idx in eachindex(ns)
        @time state = solve(NonUniformLLProblem(c=c, V=V, μ=μ, domain=(-L, L)), n_spectral=ns[idx])
        res[idx] = energy(state)
    end

    scatter(ns, abs.(res[1:end-1] .- res[end]) ./ res[end], lab="", lw=2, framestyle=:box, yscale=:log10)
end