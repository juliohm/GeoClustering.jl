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

  ctor = constructor(typeof(data))
  ctor(d, Dict(paramdim(d) => t))
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("clustering/slic.jl")
include("clustering/ghc.jl")
