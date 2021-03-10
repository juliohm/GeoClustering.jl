module GeoClustering

using Meshes
using Tables
using Distances
using Statistics
using LinearAlgebra

import Meshes: partition

include("slic.jl")

export
  SLIC,
  partition

end
