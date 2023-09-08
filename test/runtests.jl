using GeoClustering
using Meshes
using GeoTables
using GeoStatsBase
using MLJ: @load
using CategoricalArrays
using Test, Random
using Random

# list of tests
testfiles = ["mlj.jl", "slic.jl", "ghc.jl", "gsc.jl"]

@testset "GeoClustering.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end
