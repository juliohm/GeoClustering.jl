# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoClustering

using Meshes
using Tables
using TableOperations
using Distances
using Clustering
using Statistics
using LinearAlgebra
using SparseArrays

import Meshes: partition

include("clustering.jl")

export
  SLIC,
  GHC,
  partition,
  cluster

end
