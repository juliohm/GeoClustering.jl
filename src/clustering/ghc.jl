# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# auxiliary functions and variables
uniform(h; λ)      = (h ≤ λ)
triangular(h; λ)   = (h ≤ λ) * (λ - h)
epanechnikov(h; λ) = (h ≤ λ) * (λ^2 - h^2)
kernfun = Dict(
  :uniform => uniform,
  :triangular => triangular,
  :epanechnikov => epanechnikov
)

"""
    GHC(k, λ; kern=:epanechnikov, link=:ward, vars=nothing)

A method for partitioning spatial data into `k` clusters 
using Geostatistical Hierarchical Clustering (GHC).

## Parameters

* `k`    - Approximate number of clusters
* `λ`    - Approximate range of kernel function
* `kern` - Kernel function (`:uniform`, `:triangular` or `:epanechnikov`)
* `link` - Linkage function (`:single`, `:average`, `:complete`, `:ward` or `:ward_presquared`)
* `vars` - Variables (or features) to consider (default to all)

## References

* Fouedjio, F. 2016. [A hierarchical clustering method for multivariate geostatistical data]
  (https://www.sciencedirect.com/science/article/abs/pii/S2211675316300367)
"""
struct GHC <: ClusteringMethod
  k::Int
  λ::Float64
  kern::Symbol
  link::Symbol
  vars::Union{Vector{Symbol},Nothing}
end

function GHC(k, λ; kern=:epanechnikov, link=:ward, vars=nothing)
  # sanity checks
  @assert k > 0 "invalid number of clusters"
  @assert λ > 0 "invalid kernel range"
  @assert kern ∈ [:uniform,:triangular,:epanechnikov] "invalid kernel function"
  @assert link ∈ [:single,:average,:complete,:ward,:ward_presquared] "invalid linkage function"
  GHC(k, λ, kern, link, vars)
end

function partition(data, method::GHC)
  # variables used for clustering
  dvars = Tables.schema(values(data)).names
  vars  = isnothing(method.vars) ? dvars : method.vars
  @assert vars ⊆ dvars "GHC features not found in geospatial data"

  # view subset of variables
  ctor = constructor(typeof(data))
  dom  = domain(data)
  tab  = TableOperations.select(values(data), vars...)
  Ω    = ctor(dom, Dict(paramdim(dom) => tab))

  # GHC parameters
  k    = method.k
  λ    = method.λ
  kern = method.kern
  link = method.link

  # dissimilarity matrix
  D = ghc_dissimilarity_matrix(Ω, kern, λ)

  # classical hierarchical clustering
  tree = hclust(D, linkage=link)

  # cut tree to produce clusters
  labels = cutree(tree, k=k)

  # convert labels to subsets
  maxlabel = maximum(labels)
  subsets  = [Int[] for i in 1:maxlabel]
  for (i, l) in enumerate(labels)
    push!(subsets[l], i)
  end

  # return partition
  Partition(data, subsets)
end

function ghc_dissimilarity_matrix(Ω, kern, λ)
  # retrive domain/table
  𝒟 = domain(Ω)
  𝒯 = values(Ω)

  # kernel matrix
  K = ghc_kernel_matrix(kern, λ, 𝒟)

  # difference matrices
  Δ = ghc_diff_matrices(𝒯)

  # sum of cross-variograms
  Γ = ghc_variogram_sum(K, Δ)
end

function ghc_kernel_matrix(kern, λ, 𝒟)
  # kernel function
  fn    = kernfun[kern]
  Kλ(h) = fn(h, λ=λ)

  # collect coordinates
  coords = coordinates.(centroid.(𝒟))

  # lag matrix
  H = pairwise(Euclidean(), coords)

  # kernel matrix
  K = Kλ.(H)

  # return sparse version
  sparse(K)
end

function ghc_diff_matrices(𝒯)
  # covariates as columns
  covars = Tables.columntable(𝒯)

  # number of covariates
  p = length(covars)

  # one matrix per covariate pair
  Δ = Matrix{Matrix{Float64}}(undef, p, p)
  @inbounds for j in 1:p
    Zj = ghc_normalize(covars[j])
    Δj = pairwise(Euclidean(), Zj)
    for i in j+1:p
      Zi = ghc_normalize(covars[i])
      Δi = pairwise(Euclidean(), Zi)
      Δ[i,j] = Δi .* Δj
    end
    Δ[j,j] = Δj .* Δj
    for i in 1:j-1
      Δ[i,j] = Δ[j,i] # leverage the symmetry
    end
  end

  Δ
end

function ghc_variogram_sum(K, Δ)
  n = size(K, 1)
  Γ = zeros(n, n)
  for Δₒ in Δ # for each covariate pair
    # update lower triangular matrix
    @inbounds for j in 1:n
      kj = K[:,j]
      for i in j+1:n
        ki = K[:,i]
        Kij = kron(ki, kj)
        I, W = findnz(Kij)
        num = sum(W .* Δₒ[I], init=0.0)
        den = sum(W, init=0.0)
        iszero(den) || (Γ[i,j] += (1/2) * (num/den))
      end
    end
  end

  # mirror upper triangular matrix
  @inbounds for j in 1:n
    Γ[j,j] = 0.0
    for i in 1:j-1
      Γ[i,j] = Γ[j,i] # leverage the symmetry
    end
  end

  Γ
end

function ghc_normalize(x)
  μ = mean(x)
  σ = std(x, mean=μ)
  (x .- μ) ./ σ
end
