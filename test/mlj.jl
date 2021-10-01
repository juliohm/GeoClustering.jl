@testset "MLJ" begin
  Z = [ones(10,10) 2ones(10,10); 3ones(10,10) 4ones(10,10)]
  𝒮 = georef((Z=Z,))

  kmeans = @load KMeans pkg=Clustering verbosity=0
  C = cluster(𝒮, kmeans(k=4))
  @test Set(C.cluster) == Set([1,2,3,4])
end