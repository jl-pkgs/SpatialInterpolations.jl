using SpatInterp, Statistics, Test
using RTableTools
using SpatialRasterLite: st_points, st_extract

## 
@testset "interp" begin
  indir = "$(@__DIR__)/.." |> abspath
  d = fread("$indir/data/prcp_st174_shiyan.csv")

  X = [d.lon d.lat] # d.alt
  Y = d.prcp
  # Y = repeat(y, outer=(1, 24 * 30))
  b = bbox(109.5, 31.5, 112 - 0.5, 33.5)
  target = make_rast(; b, cellsize=0.01)
  # neighbor = find_neighbor(target, X; radius=100, do_angle=true)

  @time ra_idw = interp(X, Y, target; method="idw")
  @time ra_adw = interp(X, Y, target; method="adw", cdd=25, nmax=10, do_angle=true)
  @time ra_tps1 = interp_tps(X, Y, target; λ=0.01)
  @time ra_tps2 = interp_tps(X, Y, target; λ=0.1)
  # @profview ra_adw = interp(X, Y, target; wfun=weight_adw!)

  ## extract values
  points = st_points(X[:, 1], X[:, 2])
  z_idw = st_extract(ra_idw, points).value[:]
  z_adw = st_extract(ra_adw, points).value[:]
  z_tps1 = st_extract(ra_tps1, points).value[:]
  z_tps2 = st_extract(ra_tps2, points).value[:]

  @test cor(z_idw, z_adw) >= 0.88
  @test cor(z_idw, z_tps1) >= 0.88
  @test cor(z_idw, z_tps2) >= 0.88
end


## 
@testset "interp" begin
  indir = "$(@__DIR__)/.." |> abspath
  d = fread("$indir/data/prcp_st174_shiyan.csv")

  X = [d.lon d.lat] # d.alt
  Y = d.prcp

  b = bbox(109.4, 31.4, 111.6, 33.4)
  target = make_rast(; b, cellsize=1 / 120 * 2)

  ## 再加一个地理加权回归
  # neighbor = find_neighbor(target, X; radius=100, do_angle=true)
  radius = 30
  nmax = 6
  @time ra_idw = interp(X, Y, target; method="idw", radius, nmax, m=2)
  @time ra_adw = interp(X, Y, target; method="adw", radius, cdd=25, nmax, do_angle=true, m=6)
  @time ra_tps1 = interp_tps(X, Y, target; λ=0.1)
  @time ra_tps2 = interp_tps(X, Y, target; λ=0.01)
  @time ra_tps3 = interp_tps(X, Y, target; λ=0.001)

  @test true # no warning is correct
end
# @profview ra_adw = interp(X, Y, target; wfun=weight_adw!)

if false
  fig = Figure(; size=(1400, 700))
  plot_interp(fig[1, 1], ra_tps1, X, Y; title="TPS (λ = 0.1)")
  plot_interp(fig[1, 2], ra_tps2, X, Y; title="TPS (λ = 0.01)")
  plot_interp(fig[1, 3], ra_tps3, X, Y; title="TPS (λ = 0.001)")
  plot_interp(fig[2, 1], ra_idw, X, Y; title="IDW")
  plot_interp(fig[2, 2], ra_adw, X, Y; title="ADW")
  fig
end

# save("FigureS1_Prcp_2km.png", fig)
