using SpatialRasterLite, ArchGDAL, NetCDFTools
using SpatInterp
using Ipaper, DataFrames, RTableTools, JLD2, NaNStatistics
using Distances: Haversine, pairwise
using MixedGWR: BISQUARE
using ModelParams: GOF

include(joinpath(
  @__DIR__, "..", "..", "GeoWeightedRegression.jl", "example", "kfold.jl",
))
using .KfoldGWR

# # xlim = c(109.4, 111.6), ylim = c(31.2, 33.4)
# b = bbox(109.4, 31.4, 111.6, 33.4)
b = bbox(109.4, 31.2, 111.6, 33.4)
target = make_rast(; b, cellsize=1 / 120 * 2)

datadir = normpath(joinpath(@__DIR__, "..", "Project_十堰", "data"))
st = fread(joinpath(datadir, "十堰_雨量站_sp237_v20250824.csv"))
df = fread(joinpath(datadir, "十堰_prcp_sp237_mat_v20250824.csv"))

# 5.45%数据缺失
dates = DateTime.(df.time, "yyyy-mm-ddTHH:MM:SSZ")
data = Matrix(df[:, 2:end])' |> drop_missing |> collect

# 无降水不插值
_prcp = NaNStatistics.nansum(data, dims=1)[:]
inds = findall(_prcp .>= 1.0) # 所有站点总降水大于1mm日期，才进行插值
info = DataFrame(; time=dates[inds], inds, prcp=_prcp[inds])

# 交叉验证仅评估区域平均降水 >= 0.5 mm 的时次
prcp_mean = NaNStatistics.nanmean(data, dims=1)[:]
inds_rain = findall(prcp_mean .>= 0.5)

## 部分数据存在空值，如何处理
X = [st.lon st.lat] # d.alt
Y = data[:, inds]
points = map(x -> (x[1], x[2]), eachrow(X))

nlon, nlat = size(target)[1:2]
ntime = length(dates)

function spatial_weights(X_train, X_pred;
  method=:idw, radius=30, nmax=6, m=2, cdd=25,
  exclude_diagonal=false)

  ntrain, npred = size(X_train, 1), size(X_pred, 1)
  tree = SpatInterp.BallTree(collect(X_train'), Haversine(6371.0))
  W = zeros(Float64, npred, ntrain)
  nquery = min(nmax + Int(exclude_diagonal), ntrain)

  for i in axes(X_pred, 1)
    point = collect(@view X_pred[i, :])
    neighbors, distance = SpatInterp.knn(tree, point, nquery)
    order = sortperm(distance)
    neighbors, distance = neighbors[order], distance[order]

    if exclude_diagonal
      keep = neighbors .!= i
      neighbors, distance = neighbors[keep], distance[keep]
    end
    keep = distance .<= radius
    neighbors, distance = neighbors[keep], distance[keep]
    isempty(neighbors) && continue

    n = min(nmax, length(neighbors))
    neighbors, distance = neighbors[1:n], distance[1:n]
    if method === :idw
      weight = max.(distance, 1e-6) .^ (-m)
      weight ./= sum(weight)
    elseif method === :adw
      angle = angle_azimuth_sphere(point, X_train[neighbors, :])
      weight = n == 1 ? ones(1) : _weight_adw(distance, angle; cdd=Float64(cdd), m)
    else
      error("Unknown method: $method")
    end
    W[i, neighbors] .= weight
  end
  W
end

function weighted_predict(W, Y)
  valid = .!isnan.(Y)
  Y0 = ifelse.(valid, Y, 0.0)
  weight_sum = W * Float64.(valid)
  pred = W * Y0 ./ weight_sum
  pred[weight_sum .== 0] .= NaN
  pred
end

function fill_missing_idw(X, Y; radius=30, nmax=6, m=2)
  W = spatial_weights(X, X; method=:idw, radius, nmax, m,
    exclude_diagonal=true)
  Y_fill = copy(Y)
  Y_idw = weighted_predict(W, Y)
  missing = isnan.(Y_fill)
  Y_fill[missing] .= Y_idw[missing]

  # 无有效邻站时以该时次训练站平均值兜底，避免使用验证站信息
  for j in axes(Y_fill, 2)
    missing = isnan.(@view Y_fill[:, j])
    any(missing) || continue
    μ = NaNStatistics.nanmean(@view Y[:, j])
    Y_fill[missing, j] .= isnan(μ) ? 0.0 : μ
  end
  Y_fill
end

function interp_kfold(X, X_gwr, Y; k=5, seed=42, radius=30, nmax=6,
  m=2, cdd=25, λ=0.01)

  nsite = size(X, 1)
  @assert size(Y, 1) == size(X_gwr, 1) == nsite
  @assert 2 <= k <= nsite

  folds = spatial_folds(nsite; k, seed)
  methods = (:idw, :adw, :tps, :gwr_raw, :gwr_alt, :gwr_noalt, :mixed_gwr)
  pred = Dict(method => fill(NaN, size(Y)) for method in methods)
  dMat = pairwise(Haversine(6378.388), map(Tuple, eachrow(X)))

  for test in folds
    train = setdiff(axes(X, 1), test)
    X_train, X_test = X[train, :], X[test, :]
    Y_train = Y[train, :]

    W = spatial_weights(X_train, X_test; method=:idw,
      radius, nmax, m)
    pred[:idw][test, :] = weighted_predict(W, Y_train)

    W = spatial_weights(X_train, X_test; method=:adw,
      radius, nmax, m, cdd)
    pred[:adw][test, :] = weighted_predict(W, Y_train)

    Y_fill = fill_missing_idw(X_train, Y_train; radius, nmax, m)
    tps = solve_tps(X_train, Y_fill, λ)
    pred[:tps][test, :] = predict(tps, X_test; progress=false)

    p = (; np=2, adaptive=true, bw=10.0, n_max=12)
    pred[:gwr_raw][test, :] = predict_fold(
      Val(:gwr), X_gwr, dMat, train, test, Y_fill, p;
      kernel=BISQUARE, standardize=false,
    )
    pred[:gwr_noalt][test, :] = predict_fold(
      Val(:gwr), X_gwr, dMat, train, test, Y_fill, p;
      kernel=BISQUARE, standardize=true,
    )

    p = (; np=3, adaptive=true, bw=20.0, n_max=12)
    pred[:gwr_alt][test, :] = predict_fold(
      Val(:gwr), X_gwr, dMat, train, test, Y_fill, p;
      kernel=BISQUARE, standardize=true,
    )
    pred[:mixed_gwr][test, :] = predict_fold(
      Val(:mixed), X_gwr, dMat, train, test, Y_fill, p;
      kernel=BISQUARE, standardize=true,
    )
  end

  valid = .!isnan.(Y)
  for method in methods
    valid .&= .!isnan.(pred[method])
  end
  gof = Dict(method => GOF(Y[valid], pred[method][valid]) for method in methods)
  (; gof, pred, valid)
end

Y_rain = data[:, inds_rain]
dem = rast(joinpath(datadir, "dem_ShiYan_2km.tif"), FT=Float64)
alt = st_extract(dem, points).value' |> Matrix
X_gwr = [X alt]

labels = (IDW=:idw, ADW=:adw, TPS=:tps,
  GWR_raw=:gwr_raw, GWR_alt=:gwr_alt,
  GWR_noalt=:gwr_noalt, Mixed_GWR=:mixed_gwr)
cv = interp_kfold(X, X_gwr, Y_rain; k=5, seed=42)
scores = DataFrame([
  (; method=String(label), cv.gof[method]...)
  for (label, method) in pairs(labels)
])
sort!(scores, :RMSE)
@info "5-fold GOF: regional mean precipitation >= 0.5 mm" n_time=length(inds_rain) n_valid=sum(cv.valid)
display(scores)


# 以下为全场插值及结果保存，不参与交叉验证
function interpolate_all()
# TPS: 不能存在NaN，因此采用IDW插值的结果进行填补
if false
  locs = st_location(target, points)
  _inds = map(p -> LinearIndices(target.A)[p...], locs) # locs2映射到2维

  Z = nc_read("ShiYan_Prcp_Gauged237_201404-202501_2km_IDW(r=30,nmax=6).nc")
  _Z = reshape(Z, nlon * nlat, :)[_inds, :]

  data2 = deepcopy(data)
  inds_bad = isnan.(data)
  data2[inds_bad] .= _Z[inds_bad]

  jldsave("ShiYan_Pobs_interpolated_by_IDW.jld2", true; st, dates, P=data2)
end

# 缺失比例降为0.25%，余下的部分赋值为0
data2 = jldopen(joinpath(datadir, "ShiYan_Pobs_interpolated_by_IDW.jld2"))["P"]
Y2 = data2[:, inds]
Y2[isnan.(Y2)] .= 0.0 # IDW插值，有些区域无法覆盖


# TPS
for λ in [0.001, 0.01]
  @time ra_tps1 = interp_tps(X, Y, target; λ)

  R = zeros(Float32, size(target)[1:2]..., length(dates))
  R[:, :, inds] .= ra_tps1.A

  lon, lat = st_dims(target)
  dims = (; lon, lat, time=dates)
  fout = "ShiYan_Prcp_Gauged237_201404-202501_2km_TPS(λ=$λ).nc"
  ncsave(fout, true, (; units="mm h-1"); dims, P=R)
end


# IDW
if false
  radius = 30
  nmax = 6
  @time ra_idw = interp(X, Y, target; method="idw", radius, nmax, m=2)

  R = zeros(Float32, size(target)[1:2]..., length(dates))
  R[:, :, inds] .= ra_idw.A

  lon, lat = st_dims(target)
  dims = (; lon, lat, time=dates)
  fout = "ShiYan_Prcp_Gauged237_201404-202501_2km_IDW(r=30,nmax=6).nc"
  ncsave(fout, true, (; units="mm h-1"); dims, P=R)
end

# ADW
if false
  radius = 30
  nmax = 6

  for m in [2, 4, 6]
    ra = interp(X, Y, target; method="adw", radius, cdd=25, nmax, do_angle=true, m)

    R = zeros(Float32, size(target)[1:2]..., length(dates))
    R[:, :, inds] .= ra.A

    lon, lat = st_dims(target)
    dims = (; lon, lat, time=dates)
    fout = "ShiYan_Prcp_Gauged237_201404-202501_2km_ADW(r=30,nmax=6,m=$m,cdd=25).nc"
    ncsave(fout, true, (; units="mm h-1"); dims, P=R)
  end
end
end # interpolate_all
