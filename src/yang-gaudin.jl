function get_magnon_spectrum(γ, c=1.0; quadrature_rule=gausslobatto, N=100, num_points=100, kwargs...)
    rho, _, n, Q = get_ground_state(γ=γ, c=c, kwargs...)
    ε, _ = compute_dressed_energy(c, Q; kwargs...)
    xs, ws = rescale(quadrature_rule(N)..., 0., Q)
    kf = π * n

    # 2 * atan(x / (c/2))
    theta_mag(x) = 2 * atan(2 * x / c)

    # kernel is the derivative of theta
    # 1/(2pi) * 4c/(c^2+4x^2) = 2c / (pi*(c^2+4x^2))
    K_mag(x) = (2 * c) / (π * (c^2 + 4 * x^2))

    u_grid = range(0, π / 2, length=num_points)
    Λs = c .* tan.(u_grid)

    # compute P relative to the limit at Inf to fix the zero point
    # P(infinity) would be kf - dot(..., -2pi*rho) = kf + pi*n = 2kf.
    # so we define p_phys = 2kf - P_raw
    P_raw = [kf - dot(ws, (theta_mag.(xs .- Λ) .- theta_mag.(xs .+ Λ)) .* rho.(xs)) for Λ in Λs]
    E_mag = [-dot(ws, (K_mag.(xs .- Λ) .+ K_mag.(xs .+ Λ)) .* ε.(xs)) for Λ in Λs]

    # shift momentum so gapless point is at p=0
    # P_raw goes from kf (at L=0) to 2kf (at L=inf)
    # so we reverse it: p_phys goes from 0 to kf
    p_phys = abs.(P_raw[end] .- P_raw)

    return p_phys, E_mag, kf
end