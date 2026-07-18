using SpatInterp, Test
import SpatialRasterLite
using SpatialRasterLite: bbox, make_rast, st_dims
using RTableTools, Statistics


begin
  indir = "$(@__DIR__)/.." |> abspath
  d = fread("$indir/data/prcp_st174_shiyan.csv")

  X = [d.lon d.lat] # d.alt
  Y = d.prcp
  # Y = repeat(y, outer=(1, 24 * 30))
  b = bbox(109.5, 31.5, 112 - 0.5, 33.5)
  target = make_rast(; b, cellsize=0.01)
end


@testset "find_neighbor" begin
  neighbor = find_neighbor(target, X; nmax=24, radius=20, do_angle=true)
  @test mean(neighbor.count) >= 5

  neighbor = find_neighbor(target, X; nmax=24, radius=50, do_angle=true)
  @test mean(neighbor.count) >= 19
end


@testset "nearest_per_quadrant" begin
  ## nearest_per_quadrant
  neighbor = find_neighbor(target, X; nmax=24, radius=25, do_angle=true)
  neighbor4 = find_quad(neighbor)
  @test mean(neighbor4.count) >= 2.5
end

# using GLMakie, MakieLayers
# begin
#   lon, lat = st_dims(target)

#   fig = Figure(; size=(1000, 500))
#   imagesc!(fig[1, 1], lon, lat, neighbor.count)
#   imagesc!(fig[1, 2], lon, lat, neighbor4.count)
#   fig
# end
# 
