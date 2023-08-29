# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoClustering

using Meshes
using GeoStatsBase

using Tables
using Random
using TableDistances
using TableTransforms
using CategoricalArrays
using Distances
using Clustering
using Statistics
using LinearAlgebra
using SparseArrays
using ArnoldiMethod

import Meshes: partitioninds
import MLJModelInterface as MI

include("clustering.jl")

export SLIC, GHC, GSC, partition, cluster

end
