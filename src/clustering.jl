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
    cluster(data, method; vars=[all])

Cluster geospatial `data` with clustering `method`
using variables `vars` and return geospatial data
with cluster labels. Use all variables by default.

    cluster(data, model; vars=[all])

Alternatively, cluster geospatial `data` with
[`MLJ`](https://github.com/alan-turing-institute/MLJ.jl)
model.
"""
function cluster(data::Data, method; vars=nothing)
  # retrieve data
  tab = values(data)
  dom = domain(data)

  # variables used for clustering
  dvars = Tables.schema(tab).names
  vars  = isnothing(vars) ? dvars : vars
  @assert vars ⊆ dvars "variables not found in geospatial data"

  # view subset of variables
  sel = tab |> Select(vars)

  # perform clustering on selection
  _cluster(georef(sel, dom), method)
end

function _cluster(data::Data, method::ClusteringMethod)
  d = domain(data)
  p = partition(data, method)

  # assign cluster labels to samples
  labels = Vector{Int}(undef, nelements(d))
  for (l, inds) in enumerate(indices(p))
    labels[inds] .= l
  end

  # table with cluster labels
  t = (cluster=categorical(labels),)

  georef(t, d)
end

function _cluster(data::Data, model::MI.Model)
  # perform clustering
  table = values(data)
  θ, _, __ = MI.fit(model, 0, table)
  labels = MI.predict(model, θ, table)

  cluster = isprobabilistic(model) ? mode.(labels) : labels

  georef((cluster=cluster,), domain(data))
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("clustering/slic.jl")
include("clustering/ghc.jl")
