using LiebLinigerBetheAnsatz
using LinearAlgebra, Plots, LaTeXStrings

## parameters
begin
    V = x -> 0.5x^2
    L = 12.
    μ = 30.
    c = 10.
end

## default Chebyshev interpolation
begin
    @time state = solve(NonUniformLLProblem(c=2c, V=x -> 2 * V(x), μ=2μ, domain=(-L, L)), atol=1e-3)

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
    println("Chebyshev interpolation")
    println("-----------------")
    println("μ = $μ | c = $c")
    println("-----------------")
    println("Energy: $(round(e, digits=3))")
    println("Particle number: $(round(n, digits=3))")
    println("Free energy: $(round(f, digits=3))")
end
