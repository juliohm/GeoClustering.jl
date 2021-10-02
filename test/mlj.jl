@testset "MLJ" begin
  Z = [ones(10,10) 2ones(10,10); 3ones(10,10) 4ones(10,10)]
  𝒮 = georef((Z=Z,))

  # non-probabilistic model
  kmeans = @load KMeans pkg=Clustering verbosity=0
  C = cluster(𝒮, kmeans(k=4))
  @test Set(C.cluster) == Set([1,2,3,4])

  # probabilistic model
  gmm = @load GMMClusterer pkg=BetaML verbosity=0
  C = cluster(𝒮, gmm(K=4))
  # see https://github.com/alan-turing-institute/MLJ.jl/issues/846 
  @test_broken Set(C.cluster) == Set([1,2,3,4])
end