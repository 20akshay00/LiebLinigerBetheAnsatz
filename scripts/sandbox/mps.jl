using LiebLinigerBetheAnsatz
using LinearAlgebra, Plots, LaTeXStrings

begin
    μ, g = 190., 5.
    n_particles = 10
    @time state_ll = LiebLinigerBetheAnsatz.solve(FiniteLLProblem(1., 2g, μ=μ, bc=:hardwall))
    e = 0.5 * energy(state_ll)
    println(e - μ * n_particles)
    println(average_particle_density(state_ll) * 1)
end
