abstract type LLProblem end
abstract type LLState end

"""
    energy(s::LLState)
Returns the total energy E of the system. For infinite systems, returns `Inf`.
"""
function energy end

"""
    energy_density(s::LLState)
Returns the energy per unit length (E/L).
"""
function energy_density end

"""
    average_particle_density(s::LLState)
Returns the average particle density (N/L).
"""
function average_particle_density end

"""
    particle_density(s::LLState)
Returns a closure `ρ(x)` representing the spatial density profile.
"""
function particle_density end

"""
    excitation_spectrum(s::LLState; kwargs...)
Returns the elementary excitation branches: `(p_h, e_h, p_p, e_p, k_max)`.
- `p_h`, `e_h`: Momentum and Energy of the hole branch.
- `p_p`, `e_p`: Momentum and Energy of the particle branch.
- `k_max`: Maximum rapidity (analogous to kf).
"""
function excitation_spectrum end

"""
    quasimomentum_distribution(s::LLState)
Returns the distribution of quasi-momenta
"""
function quasimomentum_distribution end

"""
    fermi_quasimomentum(s::LLState)
Returns the highest quasi-momentum that is occupied.
"""
function fermi_quasimomentum end