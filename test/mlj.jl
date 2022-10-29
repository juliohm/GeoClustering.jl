@testset "MLJ" begin
  Ω = georef((Z=[ones(10,10) 2ones(10,10); 3ones(10,10) 4ones(10,10)],))

  # non-probabilistic model
  m = @load KMeans pkg=Clustering verbosity=0
  C = cluster(Ω, m(k=4))
  @test Set(C.cluster) == Set(categorical([1,2,3,4]))

  # probabilistic model
  m = @load GaussianMixtureClusterer pkg=BetaML verbosity=0
  C = cluster(Ω, m(n_classes=4))
  @test Set(C.cluster) == Set(categorical([1,2,3,4]))
end