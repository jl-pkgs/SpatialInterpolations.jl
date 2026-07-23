using GLMakie, MakieLayers
include("main_vis.jl")

radius = 30
nmax = 6
@time ra_idw = interp(X, Y, target; method="idw", radius, nmax, m=2)
@time ra_adw = interp(X, Y, target; method="adw", radius, cdd=25, nmax, do_angle=true, m=6)
@time ra_tps1 = interp_tps(X, Y, target; λ=0.1)
@time ra_tps2 = interp_tps(X, Y, target; λ=0.01)
@time ra_tps3 = interp_tps(X, Y, target; λ=0.001)
# @profview ra_adw = interp(X, Y, target; wfun=weight_adw!)

fig = Figure(; size=(1400, 700))
plot_interp(fig[1, 1], ra_idw, X, Y; title="IDW")
plot_interp(fig[1, 2], ra_adw, X, Y; title="ADW")
# plot_interp(fig[2, 1], ra_tps1, X, Y; title="TPS (λ = 0.1)")
# plot_interp(fig[2, 2], ra_tps2, X, Y; title="TPS (λ = 0.01)")
# plot_interp(fig[2, 3], ra_tps3, X, Y; title="TPS (λ = 0.001)")
fig
