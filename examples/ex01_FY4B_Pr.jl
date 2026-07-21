using NetCDFTools: nc_read, st_dims
using CairoMakie, MakieLayers
using SpatInterp, JLD2
using SpatialRasterLite: SpatRaster, bbox, make_rast, bbox2dims, rast

## 
import MakieLayers: imagesc, imagesc!
function imagesc!(handle, ra::SpatRaster, args...; kw...)
  lon, lat = st_dims(ra)
  imagesc!(handle, lon, lat, ra.A, args...; kw...)
end

f = "data/FY4B_20250714014500_20250714015959.nc"
# lon, lat = st_dims(f)

I, J = nc_read(f, "i"), nc_read(f, "j")
lon = nc_read(f, "lon")
lat = nc_read(f, "lat")
Prcp = nc_read(f, "P")

# 采用bilinear_irregular，对数据进行插值
inds = findall(.!isnan.(lon[:]))
X = [lon[inds] lat[inds]] # 76.6%, lon, lat
Z = repeat(Prcp[inds], outer=(1, 10))

function make_rast2(; b::bbox, cellsize, FT=Float64)
  lon, lat = bbox2dims(b; cellsize)
  nlon, nlat = length(lon), length(lat)
  rast(zeros(FT, nlon, nlat), b)
end

b = bbox(70., 15., 140., 55.)
target = make_rast2(; b, cellsize=1 / 120 * 4, FT=Float32)

# radius = 4km * 4
fast = true
cellsize = 1 / 120 * 4 # ~4km,
radius = fast ? cellsize * 10 : 4 * 5

@time neighbor = find_neighbor(target, X; nmax=36, fast, radius, do_angle=true)
neighbor4 = find_quad(neighbor)
@time res = interp_weight(neighbor4, Z)

# sum(isnan.(res.A)) / prod(size(res)) # 27.7%的NaN

begin
  fig = Figure(; size=(1500, 600))
  imagesc!(fig[1, 1], I, J, Prcp; title="P")
  imagesc!(fig[1, 2], res[:, :, 1:4])
  fig
end

## 选择其中一个进行查看
sum(isnan.(neighbor4.weight)) / prod(size(res)) # 27.7%的NaN

begin
  (; count, weight, dist) = neighbor4
  fig = Figure(; size=(1500, 600))
  # imagesc!(fig[1, 1], I, J, count; title="count")
  # imagesc!(fig, I, J, weight; title="count")
  # imagesc!(fig, I, J, dist; title="dist")
  # imagesc!(fig, I, J, neighbor4.angle; title="angle")
  # imagesc!(fig[1, 2], I, J, lat; title="lat")
  # imagesc!(fig[1, 3], I, J, P; title="P")
  fig
end

begin
  fig = Figure(; size=(1500, 600))
  imagesc!(fig[1, 1], I, J, lon; title="lon")
  imagesc!(fig[1, 2], I, J, lat; title="lat")
  imagesc!(fig[1, 3], I, J, P; title="P")
  fig
end
