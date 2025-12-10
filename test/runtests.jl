using DrWatson, Test
@quickactivate "LiebLinigerBetheAnsatz"
using LiebLinigerBetheAnsatz, LinearAlgebra
# include(srcdir("file.jl"))

println("Starting tests")

@testset "Fredholm solvers" begin
    xlim = [0, 1]
    xs = range(xlim..., 100)
    f_exact(x) = exp(x)

    f1, _, _ = solve(
        QuadratureSolver(gausslegendre(30)),
        (x, y) -> -min(x, y),
        x -> x * ℯ + 1,
        xlim...
    )

    f2, _, _ = solve(
        ModifiedQuadratureSolver(gausslegendre(30)),
        (x, y) -> -min(x, y),
        x -> x * ℯ + 1,
        x -> -(x - x^2 / 2),
        xlim...
    )

    @test sum(abs.(f1.(xs) .- f_exact.(xs))) / length(xs) < 1e-3
    @test sum(abs.(f2.(xs) .- f_exact.(xs))) / length(xs) < 1e-5

    xlim = [-1., 1.]
    xs = range(xlim..., 100)
    function f_exact(x)
        c2 = (ℯ^4 + 6 * ℯ^2 + 1) / (8 * (ℯ^2 + 1))
        c1 = c2 + 1 / (1 + ℯ^2)
        return (c1 + 0.5 * x) * exp(x) + c2 * exp(-x)
    end

    f1, _, _ = solve(
        QuadratureSolver(gausslegendre(30)),
        (x, y) -> 0.5 * abs(x .- y),
        x -> exp(x),
        xlim...
    )

    f2, _, _ = solve(
        ModifiedQuadratureSolver(gausslegendre(30)),
        (x, y) -> 0.5 * abs(x .- y),
        x -> exp(x),
        (x) -> 0.5 * (1 + x^2),
        xlim...
    )

    @test sum(abs.(f1.(xs) .- f_exact.(xs))) / length(xs) < 1e-1
    @test sum(abs.(f2.(xs) .- f_exact.(xs))) / length(xs) < 1e-4
end

@testset "Lieb-Liniger ground state" begin
    # Tonks-Girardeau limit
    rho, e, n, Q = get_ground_state(γ=1e8, N=100)
    @test norm(e / n^3 - π^2 / 3) < 1e-6
    @test norm(n - Q / pi) < 1e-10
    @test norm(rho(0) - 1 / 2π) < 1e-8
end
