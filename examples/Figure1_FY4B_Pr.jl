using NetCDFTools: nc_read, st_dims
using SpatInterp, JLD2
using SpatialRasterLite: SpatRaster, bbox, make_rast, bbox2dims, rast
using GLMakie, MakieLayers, Shapefile
using Shapefile: Table

shp = Table("Z:/Global/GlobalWaterBalance/GlobalWB/data/shp/GlobalLand.shp")
poly_china = Table("D:/Documents/GitHub/nmc_met_graphics/nmc_met_graphics/resources/maps/bou2_4p.shp")
arc_china = Table("D:/Documents/GitHub/nmc_met_graphics/nmc_met_graphics/resources/maps/bou1_4l.shp")

import MakieLayers: imagesc, imagesc!
function imagesc!(handle, ra::SpatRaster, args...; kw...)
  lon, lat = st_dims(ra)
  imagesc!(handle, lon, lat, ra.A, args...; kw...)
end

function make_rast2(; b::bbox, cellsize, FT=Float64)
  lon, lat = bbox2dims(b; cellsize)
  nlon, nlat = length(lon), length(lat)
  rast(zeros(FT, nlon, nlat), b)
end



begin
  f = "data/FY4B_20250714014500_20250714015959.nc"

  I, J = nc_read(f, "i"), nc_read(f, "j")
  lon, lat = st_dims(f)
  Prcp = nc_read(f, "P")

  b = bbox(70., 15., 140., 55.)
  inds = @.((b.xmin <= lon <= b.xmax) & (b.ymin <= lat <= b.ymax)) |> findall
  # inds = findall(.!isnan.(lon[:]))

  # 采用bilinear_irregular，对数据进行插值
  X = [lon[inds] lat[inds]] # 76.6%, lon, lat
  Z = Prcp[inds]
  # Z = repeat(Prcp[inds], outer=(1, 10))
  target = make_rast2(; b, cellsize=1 / 120 * 4, FT=Float32)
end


begin
  # radius = 4km * 4
  fast = true
  cellsize = 1 / 120 * 4 # ~4km,
  radius = fast ? cellsize * 10 : 4 * 5

  @time neighbor = find_neighbor(target, X; nmax=36, fast, radius, do_angle=true)
  neighbor4 = find_quad(neighbor)
  @time res = interp_weight(neighbor4, Z)
end

Z = Prcp[inds]
@time res = interp_weight(neighbor4, Z)
# Z = repeat(Prcp[inds], outer=(1, 1000)) # test speed
# @time res = interp_weight(neighbor4, Z)

begin
  fig = Figure(; size=(1500, 550))
  imagesc!(fig[1, 1], I, J, Prcp; title="RAW")

  ax, plt = imagesc!(fig[1, 2], res; title="Irregular Bilinear")
  lines!(ax, arc_china.geometry, linewidth=0.5, alpha=1, color=:black)
  poly!(ax, poly_china.geometry, color=nan_color, strokewidth=0.2, strokecolor=:black)
  poly!(ax, shp.geometry, color=nan_color, strokewidth=0.5, strokecolor=:black)

  set_lims!([ax], (70, 140), (15, 55))
  fig

  save("Figure1_FY4B_Pr.png", fig)
end
