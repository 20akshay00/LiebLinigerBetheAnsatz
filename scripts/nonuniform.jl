using LiebLinigerBetheAnsatz
using LinearAlgebra, Plots, LaTeXStrings

## parameters
begin
    V = x -> 0.5x^2
    L = 5.
    μ = 5.
    c = 10.
end

## default Chebyshev interpolation
begin
    @time state = solve(NonUniformLLProblem(c=c, V=V, μ=μ, domain=(-L, L)), atol=1e-3)

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

## custom interpolant
begin
    # linear interpolation
    function lerp(mu, x_grid, y_vals)
        idx = searchsortedfirst(x_grid, mu)
        idx <= 1 && return 0.0
        idx > length(x_grid) && return y_vals[end]

        x0, x1 = x_grid[idx-1], x_grid[idx]
        y0, y1 = y_vals[idx-1], y_vals[idx]

        return y0 + (y1 - y0) * (mu - x0) / (x1 - x0)
    end

    function construct_lerp_eos(c::Float64, μ_max::Float64; L=200)
        μs = range(1e-8, μ_max, length=L)
        states = [solve(InfiniteLLProblem(c=c, μ=μ)) for μ in μs]
        ρs = average_particle_density.(states)
        es = energy_density.(states)
        return (μ -> lerp(μ, μs, ρs), μ -> lerp(μ, μs, es))
    end
end

begin
    ρ_eos, e_eos = construct_lerp_eos(c, μ, L=200)
    state1 = solve(NonUniformLLProblem(c=c, V=V, μ=μ, ρ_eos=ρ_eos, e_eos=e_eos, domain=(-L, L)))
    dens1 = particle_density(state1)
    xs1 = range(domain(state1)..., 100)
    plot(dens1, xs1, lab="E = $(round(energy(state1), digits=3))", lw=2)
    plot!(xlabel="Position", ylabel="Single particle density", framestyle=:box)
end

begin
    e1 = energy(state1) * 0.5
    n1 = sum(dens1.(xs) * step(xs1))
    f1 = e1 - μ * n1

    println("-----------------")
    println("Linear interpolation")
    println("-----------------")
    println("μ = $μ | c = $c")
    println("-----------------")
    println("Energy: $(round(e1, digits=3))")
    println("Particle number: $(round(n1, digits=3))")
    println("Free energy: $(round(f1, digits=3))")
end
