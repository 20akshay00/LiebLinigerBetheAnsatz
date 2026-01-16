"""
    get_magnon_spectrum(γ, c=1.0; quadrature_rule=gausslobatto, N=100, num_points=100, kwargs...)

Calculate the magnon excitation spectrum (dark solitons) for the Lieb-Liniger model using 
the dressed energy and phase shift formalism.

# Arguments
- `γ`: Dimensionless interaction strength (c/n).
- `c`: Coupling constant.
- `quadrature_rule`: Integration rule for the solvers.
- `N`: Number of quadrature points.
- `num_points`: Number of points used to sample the spectrum.
- `kwargs`: Keyword arguments passed to ground state and dressed energy solvers.

# Returns
- `p_phys`: Physical momentum of the magnon branch.
- `e_mag`: Excitation energy of the magnon branch.
- `kf`: The Fermi momentum π * n.
"""
function get_magnon_spectrum(γ, c=1.0; quadrature_rule=gausslobatto, N=100, num_points=100, kwargs...)
    rho, _, n, Q = get_ground_state(γ=γ, c=c, kwargs...)
    ε, _ = compute_dressed_energy(c, Q; kwargs...)
    xs, ws = rescale(quadrature_rule(N)..., 0., Q)
    kf = π * n

    # 2 * atan(x / (c/2))
    θ(x) = 2 * atan(2 * x / c)

    # kernel is the derivative of theta
    # 1/(2pi) * 4c/(c^2+4x^2) = 2c / (pi*(c^2+4x^2))
    K(x) = (2 * c) / (π * (c^2 + 4 * x^2))

    Λs = c .* tan.(range(0, π / 2, length=num_points))

    # compute P relative to the limit at Inf to fix the zero point
    # P(infinity) would be kf - dot(..., -2pi*rho) = kf + pi*n = 2kf.
    # so we define p_phys = 2kf - P_raw
    p_raw = [kf - dot(ws, (θ.(xs .- Λ) .- θ.(xs .+ Λ)) .* rho.(xs)) for Λ in Λs]
    e_mag = [-dot(ws, (K.(xs .- Λ) .+ K.(xs .+ Λ)) .* ε.(xs)) for Λ in Λs]

    # shift momentum so gapless point is at p=0
    # P_raw goes from kf (at L=0) to 2kf (at L=inf)
    # so we reverse it: p_phys goes from 0 to kf
    p_phys = abs.(p_raw[end] .- p_raw)

    return p_phys, e_mag, kf
end