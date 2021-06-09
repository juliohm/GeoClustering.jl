# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoClustering

using Meshes
using Tables
using TableOperations
using Distances
using Statistics
using LinearAlgebra

import Meshes: partition

include("clustering.jl")

export
  SLIC,
  partition,
  cluster

end
