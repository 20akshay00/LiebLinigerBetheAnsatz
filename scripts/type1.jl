using LiebLinigerBetheAnsatz, Plots, LaTeXStrings

γ_vals = [0.5, 5.0, 50.0, 500.0]
N_pts = 100
# Fix limits so steepness is visually comparable
xlims_global = (0, 2π) 
ylims_global = (0, 35) # High enough to see curvature, low enough for TG detail

plots = [begin
    rho, e, n, Q = get_ground_state(γ=γ, N=N_pts)
    p_h, E_h, p_p, E_p = get_particle_hole_spectrum(γ, N=N_pts, num_points=40)
    
    # Reference Bogoliubov
    E_bdg(p) = sqrt(p^4 + 4*γ*n^2*p^2)
    p_grid = range(0, 2π*n, length=100)
    
    plt = plot(title=L"\gamma = %$γ", framestyle=:box, grid=false)
    
    # Unified Scaling: Energy / n^2 and Momentum / n
    plot!(p_grid ./ n, E_bdg.(p_grid) ./ n^2, ls=:dash, c=:black, alpha=0.6, lab=(γ==γ_vals[1] ? "BdG" : ""))
    scatter!(p_p ./ n, E_p ./ n^2, label=(γ==γ_vals[1] ? "Type I" : ""), ms=2.5, c=1, markerstrokewidth=0)
    scatter!(p_h ./ n, E_h ./ n^2, label=(γ==γ_vals[1] ? "Type II" : ""), ms=2.5, c=2, markerstrokewidth=0)
    
    xlims!(xlims_global...)
    ylims!(ylims_global...)
    
    xlabel!(L"p/n")
    if i == 1 ylabel!(L"\epsilon/n^2") end # Only label y-axis on first plot
    plt
end for (i, γ) in enumerate(γ_vals)]

plot(plots..., layout=(1, 4), size=(1200, 350), left_margin=5Plots.mm)