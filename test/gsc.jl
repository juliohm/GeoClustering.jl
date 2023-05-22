@testset "GSC" begin
  Ω = georef((Z=[10sin(i / 10) + j for i in 1:100, j in 1:100],))

  Random.seed!(2022)
  C = cluster(Ω, GSC(50, 2.0))
  @test Set(C.cluster) == Set(1:50)

  if visualtests
    @test_reference "data/gsc.png" plot(C)
  end
end
