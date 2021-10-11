# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SLIC(k, m; tol=1e-4, maxiter=10)

A method for clustering geospatial data into approximately `k`
clusters using Simple Linear Iterative Clustering (SLIC).
The method produces clusters of samples that are spatially
connected based on a distance `dₛ` and that, at the same
time, are similar in terms of `vars` with distance `dᵥ`.
The tradeoff is controlled with a hyperparameter parameter
`m` in an additive model `dₜ = √(dᵥ² + m²(dₛ/s)²)`.

## Parameters

* `k`       - Approximate number of clusters
* `m`       - Hyperparameter of SLIC model
* `tol`     - Tolerance of k-means algorithm (default to `1e-4`)
* `maxiter` - Maximum number of iterations (default to `10`)

## References

* Achanta et al. 2011. [SLIC superpixels compared to state-of-the-art
  superpixel methods](https://ieeexplore.ieee.org/document/6205760)
"""
struct SLIC <: ClusteringMethod
  k::Int
  m::Float64
  tol::Float64
  maxiter::Int
end

function SLIC(k::Int, m::Real; tol=1e-4, maxiter=10)
  @assert tol > 0 "invalid tolerance"
  @assert maxiter > 0 "invalid number of iterations"
  SLIC(k, m, tol, maxiter)
end

function partition(data, method::SLIC)
  # normalize atributes
  𝒯 = TableDistances.normalize(values(data))
  Ω = georef(first(𝒯), domain(data))

  # SLIC hyperparameter
  m = method.m

  # initial spacing of clusters
  s = slic_spacing(Ω, method)

  # initialize cluster centers
  c = slic_initialization(Ω, s)

  # ball neighborhood search
  searcher = BallSearch(Ω, NormBall(s))

  # pre-allocate memory for label and distance
  l = fill(0, nelements(Ω))
  d = fill(Inf, nelements(Ω))

  # performance parameters
  tol     = method.tol
  maxiter = method.maxiter

  # k-means algorithm
  err, iter = Inf, 0
  while err > tol && iter < maxiter
    o = copy(c)

    slic_assignment!(Ω, searcher, m, s, c, l, d)
    slic_update!(Ω, c, l)

    err = norm(c - o) / norm(o)
    iter += 1
  end

  subsets = [findall(isequal(k), l) for k in 1:length(c)]

  Partition(data, subsets)
end

function slic_spacing(data, method)
  V = measure(boundingbox(data))
  d = embeddim(data)
  k = method.k
  (V/k) ^ (1/d)
end

function slic_initialization(data, s)
  # efficient neighbor search
  searcher = KNearestSearch(data, 1)

  # bounding box properties
  bbox = boundingbox(data)
  lo, up = coordinates.(extrema(bbox))

  # cluster centers
  clusters = Vector{Int}()
  neighbor = Vector{Int}(undef, 1)
  ranges = [(l+s/2):s:u for (l, u) in zip(lo, up)]
  for x in Iterators.product(ranges...)
    search!(neighbor, Point(x), searcher)
    push!(clusters, neighbor[1])
  end

  unique(clusters)
end

function slic_assignment!(data, searcher, m, s, c, l, d)
  for (k, cₖ) in enumerate(c)
    pₖ = centroid(data, cₖ)
    inds = search(pₖ, searcher)

    # distance between points
    X  = (coordinates(centroid(data, ind)) for ind in inds)
    xₖ = [coordinates(pₖ)]
    dₛ = pairwise(Euclidean(), X, xₖ)

    # distance between variables
    𝒮ᵢ = view(data, inds)
    𝒮ₖ = view(data, [cₖ])
    V  = values(𝒮ᵢ)
    vₖ = values(𝒮ₖ)
    dᵥ = pairwise(TableDistance(normalize=false), V, vₖ)

    # total distance
    dₜ = @. √(dᵥ^2 + m^2 * (dₛ/s)^2)

    @inbounds for (i, ind) in enumerate(inds)
      if dₜ[i] < d[ind]
        d[ind] = dₜ[i]
        l[ind] = k
      end
    end
  end
end

function slic_update!(data, c, l)
  for k in 1:length(c)
    inds = findall(isequal(k), l)
    X  = (coordinates(centroid(data, ind)) for ind in inds)
    μ  = [mean(X)]
    dₛ = pairwise(Euclidean(), X, μ)
    @inbounds c[k] = inds[argmin(vec(dₛ))]
  end
end
