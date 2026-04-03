using DrWatson, Test
@quickactivate "LiebLinigerBetheAnsatz"
using LiebLinigerBetheAnsatz, LinearAlgebra

# helper to extract standard observables for tests
function get_observables(s)
    return (
        quasimomentum_distribution(s),
        energy_density(s),
        average_particle_density(s),
        fermi_quasimomentum(s)
    )
end

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

@testset "Lieb-Liniger solver suite" begin
    @testset "Infinite system limits" begin
        # TG limit
        # exact: e/n^3 = π²/3, Q = π*n, ρ(k) = 1/2π
        target_n = 1.0
        prob_tg = InfiniteLLProblem(c=1e8, γ=1e8) # high gamma
        state_tg = solve(prob_tg)
        rho, e, n, Q = get_observables(state_tg)

        @test e / n^3 ≈ π^2 / 3 atol = 1e-6
        @test Q ≈ π * n atol = 1e-6
        @test rho(0) ≈ 1 / (2π) atol = 1e-6

        # BEC limit
        # exact: e -> 0, Q -> π*n (but rho is peaked at 0)
        prob_weak = InfiniteLLProblem(c=1e-5, μ=0.1)
        state_weak = solve(prob_weak)
        _, e_w, n_w, _ = get_observables(state_weak)

        @test e_w >= 0.0
        #@test e_w < 1e-3 # should be small, not sure how small
    end

    @testset "Finite system PBC" begin
        L = 10.0
        for N in [2, 10, 20]
            # TG limit
            prob_tg = FiniteLLProblem(L, 1e8; N=N, bc=:periodic)
            state_tg = solve(prob_tg)
            k_solved = quasimomentum_distribution(state_tg)

            ns = collect(1:N) .- (N + 1) / 2
            k_expected = (2π / L) .* ns
            @test k_solved ≈ k_expected atol = 1e-5

            # BEC limit
            prob_bec = FiniteLLProblem(L, 1e-10; N=N, bc=:periodic)
            state_bec = solve(prob_bec)
            @test all(abs.(quasimomentum_distribution(state_bec)) .< 1e-4)
        end
    end

    @testset "Finite system OBC" begin
        L = 10.0
        for N in [2, 10, 20]
            # TG limit
            # k_j = j * π / L
            prob_tg = FiniteLLProblem(L, 1e8; N=N, bc=:hardwall)
            state_tg = solve(prob_tg)
            k_solved = sort(quasimomentum_distribution(state_tg))

            k_expected = [j * π / L for j in 1:N]
            @test k_solved ≈ k_expected atol = 1e-4

            # BEC limit
            # All bosons collapse into the n=1 mode: k = π/L
            # E_total = N * (π/L)²
            prob_bec = FiniteLLProblem(L, 0.0; N=N, bc=:hardwall)
            state_bec = solve(prob_bec)
            k_bec = quasimomentum_distribution(state_bec)
            e_bec = energy_density(state_bec) * L

            @test all(abs.(k_bec .- π / L) .< 1e-4)
            @test e_bec ≈ N * (π / L)^2 atol = 1e-5
        end
    end

    @testset "Thermodynamic Convergence" begin
        # test if finite system energy density approaches TDL energy density
        c = 1.0
        rho_target = 1.0

        # TDL value
        prob_inf = InfiniteLLProblem(c=c, γ=c / rho_target)
        state_inf = solve(prob_inf)
        e_tdl = energy_density(state_inf)

        # Finite system (increasing N)
        N = 500
        L = N / rho_target
        prob_fin = FiniteLLProblem(L, c; N=N, bc=:periodic)
        state_fin = solve(prob_fin)
        e_finite = energy_density(state_fin)

        println(e_finite - e_tdl)
        @test e_finite ≈ e_tdl rtol = 1e-5
    end

    @testset "Numerical Stability & Convexity" begin
        # test large N convergence with damping
        N = 1000
        L = 1000.0
        c = 0.5
        prob = FiniteLLProblem(L, c; N=N, bc=:periodic)

        # should not throw any errors
        state = solve(prob; relaxation=0.5, maxiter=100)
        @test length(quasimomentum_distribution(state)) == N
    end
end