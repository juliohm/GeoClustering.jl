@testset "GSC" begin
  Random.seed!(2022)
  Ω = georef((Z=[10sin(i / 10) + j for i in 1:100, j in 1:100],))
  C = cluster(Ω, GSC(50, 2.0))
  @test Set(C.cluster) == Set(1:50)
end
