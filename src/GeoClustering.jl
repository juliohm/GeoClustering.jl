# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoClustering

using Meshes
using GeoTables
using GeoStatsBase

using Random
using Tables
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
