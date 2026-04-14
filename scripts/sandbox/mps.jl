using LiebLinigerBetheAnsatz
using LinearAlgebra, Plots, LaTeXStrings

begin
    μ, g = 160., 10.
    xmax = 0.5
    @time state_ll = LiebLinigerBetheAnsatz.solve(FiniteLLProblem(2xmax, 2g, μ=2μ, bc=:hardwall))
    n_particles = state_ll.N
    println(n_particles, " particles.")
    e = 0.5 * energy(state_ll)
    println(e - μ * n_particles)
end
