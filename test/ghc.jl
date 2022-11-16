@testset "GHC" begin
  Z = [ones(10,10) 2ones(10,10); 3ones(10,10) 4ones(10,10)]
  𝒮 = georef((Z=Z,))
  p = partition(𝒮, GHC(4, 1.0))
  @test length(p) == 4
  d1, d2, d3, d4 = domain.(p)
  @test all(nelements.([d1,d2,d3,d4]) .== 100)
  @test mean(coordinates.(centroid.(d1))) == [5.0,5.0]
  @test mean(coordinates.(centroid.(d2))) == [15.0,5.0]
  @test mean(coordinates.(centroid.(d3))) == [5.0,15.0]
  @test mean(coordinates.(centroid.(d4))) == [15.0,15.0]

  C = cluster(𝒮, GHC(4, 1.0))
  @test C.cluster == vec(Z')

  𝒮′ = georef(values(𝒮), centroid.(domain(𝒮)))
  C′ = cluster(𝒮′, GHC(4, 1.0))
  @test C.cluster == C′.cluster

  𝒮 = georef((z=[√(i^2+j^2) for i in 1:50, j in 1:50],))
  p = partition(𝒮, GHC(50, 1.0))
  @test length(p) == 50

  if visualtests
    @test_reference "data/ghc.png" plot(p)
  end
end
