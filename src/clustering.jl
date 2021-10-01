# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ClusteringMethod

A method for clustering geospatial data.

Unlike a generic [`PartitionMethod`](@ref), a
clustering method also uses variables in the
data besides the underlying geospatial domain.
"""
abstract type ClusteringMethod <: Meshes.PartitionMethod end

"""
    cluster(data, method)

Cluster geospatial `data` with clustering `method`
and return geospatial data with cluster labels.
"""
function cluster(data::Data, method::ClusteringMethod)
  d = domain(data)
  p = partition(data, method)

  # assign cluster labels to samples
  labels = Vector{Int}(undef, nelements(d))
  for (l, inds) in enumerate(indices(p))
    labels[inds] .= l
  end

  # table with cluster labels
  t = (cluster = labels,)

  georef(t, d)
end

"""
    cluster(data, model)

Cluster geospatial `data` with model implementing the
[`MLJ`](https://github.com/alan-turing-institute/MLJ.jl)
interface.
"""
function cluster(data::Data, model::MI.Model)
  # perform clustering
  table = values(data)
  θ, _, __ = MI.fit(model, 0, table)
  labels = MI.predict(model, θ, table)

  # table with cluster labels
  t = (cluster = unwrap.(labels),)

  # underlying geospatial domain
  d = domain(data)

  georef(t, d)
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("clustering/slic.jl")
include("clustering/ghc.jl")
