using ModelParams: GOF
using Random: MersenneTwister, randperm
using RTableTools
using SpatialRasterLite: st_extract
using SpatInterp, Statistics, Test

function interp_kfold(X, y, target; k=5, seed=42)
  n = length(y)
  @assert size(X, 1) == n
  @assert 2 <= k <= n

  order = randperm(MersenneTwister(seed), n)
  folds = [order[i:k:end] for i in 1:k]
  methods = (:idw, :adw, :tps)
  pred = Dict(method => fill(NaN, n) for method in methods)
  points = st_points(X[:, 1], X[:, 2])

  for test in folds
    train = setdiff(eachindex(y), test)
    X_train, y_train = X[train, :], y[train]

    ra = interp(X_train, y_train, target;
      method="idw", radius=200, nmax=10, m=2)
    pred[:idw][test] = st_extract(ra, points[test]).value[:]

    ra = interp(X_train, y_train, target;
      method="adw", radius=200, nmax=10,
      do_angle=true, cdd=25, m=2)
    pred[:adw][test] = st_extract(ra, points[test]).value[:]

    tps = solve_tps(X_train, y_train, 0.01)
    pred[:tps][test] = predict(tps, X[test, :]; progress=false)[:]
  end

  Dict(method => GOF(y, pred[method]) for method in methods)
end

@testset "5-fold interpolation" begin
  indir = joinpath(@__DIR__, "..")
  d = fread(joinpath(indir, "data", "prcp_st174_shiyan.csv"))
  X = [d.lon d.lat]
  y = Float64.(d.prcp)

  cellsize = 0.01
  b = bbox(
    minimum(X[:, 1]) - cellsize,
    minimum(X[:, 2]) - cellsize,
    maximum(X[:, 1]) + cellsize,
    maximum(X[:, 2]) + cellsize,
  )
  target = make_rast(; b, cellsize)
  scores = interp_kfold(X, y, target; k=5)

  for method in (:idw, :adw, :tps)
    gof = scores[method]
    @info "5-fold GOF" method gof
    @testset "$method" begin
      @test gof.n_valid == length(y)
      @test gof.NSE > 0.4
      @test gof.R2 > 0.5
      @test gof.KGE > 0.6
      @test gof.RMSE < std(y)
    end
  end
end
