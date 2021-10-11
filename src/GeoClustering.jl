# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoClustering

using Meshes
using GeoStatsBase

using Tables
using TableDistances
using TableOperations
using CategoricalArrays
using Distances
using Clustering
using Statistics
using LinearAlgebra
using SparseArrays

import Meshes: partition
import MLJModelInterface as MI

include("clustering.jl")

export
  SLIC,
  GHC,
  partition,
  cluster

end
