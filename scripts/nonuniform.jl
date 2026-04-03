using LiebLinigerBetheAnsatz
using LinearAlgebra, Plots, LaTeXStrings

begin
    V = x -> 0.5x^2
    L = 7.
    μ = 10.
    c = 1e-1
    state = solve(NonUniformLLProblem(c=c, V=V, μ=μ, domain=(-L, L)))
    dens = particle_density(state)
    xs = range(domain(state)..., 100)
    plot(dens, xs, lab="E = $(round(energy(state), digits=3))", lw=2)
    plot!(xlabel="Position", ylabel="Single particle density", framestyle=:box)
end

begin
    e = energy(state)
    n = sum(dens.(xs) * step(xs))
    f = e - μ * n

    println("-----------------")
    println("μ = $μ | c = $c")
    println("-----------------")
    println("Energy: $(round(e, digits=3))")
    println("Particle number: $(round(n, digits=3))")
    println("Free energy: $(round(f, digits=3))")
end