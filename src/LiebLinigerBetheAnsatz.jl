module LiebLinigerBetheAnsatz

using Reexport
@reexport using FastGaussQuadrature
using Roots, LinearAlgebra, LaTeXStrings

include("fredholm.jl")
include("lieb-liniger.jl")
include("utils.jl")
include("yang-gaudin.jl")

export QuadratureSolver, ModifiedQuadratureSolver, solve
export get_ground_state, get_particle_hole_spectrum, get_magnon_spectrum
export pitick, reimann_quadrature, midpoint_quadrature, trapezoidal_quadrature, simpson_quadrature
export compute_chemical_potential, compute_dressed_energy

end