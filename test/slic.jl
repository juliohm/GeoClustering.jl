@testset "SLIC" begin
  Z = [ones(10,10) 2ones(10,10); 3ones(10,10) 4ones(10,10)]
  𝒮 = georef((Z=Z,))
  p = partition(𝒮, SLIC(4, 1.0))
  @test length(p) == 4
  @test all(nelements.(p) .== 100)
  p1, p2, p3, p4 = p
  @test mean(coordinates(centroid(p1, ind)) for ind in 1:nelements(p1)) == [5.0,5.0]
  @test mean(coordinates(centroid(p2, ind)) for ind in 1:nelements(p2)) == [15.0,5.0]
  @test mean(coordinates(centroid(p3, ind)) for ind in 1:nelements(p3)) == [5.0,15.0]
  @test mean(coordinates(centroid(p4, ind)) for ind in 1:nelements(p4)) == [15.0,15.0]

  C = cluster(𝒮, SLIC(4, 1.0))
  @test C.cluster == vec(Z')

  𝒮 = georef((z=[√(i^2+j^2) for i in 1:100, j in 1:100],))
  p = partition(𝒮, SLIC(50, 0.001))
  @test 50 ≤ length(p) ≤ 60

  if visualtests
    @test_reference "data/slic.png" plot(p)
  end

  # test SLIC with heterogeneous data
  Z = (a=rand(10), b=1:10, x=rand(10), y=rand(10))
  𝒮 = georef(Z, (:x, :y))
  C = cluster(𝒮, SLIC(2, 1.0))
  @test domain(C) == domain(𝒮)
  @test Set(C.cluster) ⊆ Set(1:2)

  # test SLIC for orphaned points
  a = [0,0,0,0,0,0,0,0,0,0]
  x = [0.4993029939801461, 0.14954882636793432, 0.23118957975519616, 0.6816610871344635, 0.6665309965318731, 0.691522274292691, 0.012495903053589608, 0.9831177095525963, 0.4445263730141056, 0.2175871587746574]
  y = [0.32721108209880256, 0.11427387079564899, 0.826401075107011, 0.6164294766961782, 0.6562529361193601, 0.43388375115444644, 0.7624847842129086, 0.1516623758764959, 0.07641616063237144, 0.8669098569279463]
  Z = (a=a, x=x, y=y)
  𝒮 = georef(Z, (:x, :y))
  C = cluster(𝒮, SLIC(2, 1.0))
  @test Set(C.cluster) ⊆ Set(1:2)

  # test SLIC with weights in attribute columns
  z1 = [√((i-0)^2+(j-0)^2) for i in 1:100, j in 1:100]
  z2 = [√((i-100)^2+(j-100)^2) for i in 1:100, j in 1:100]
  𝒮 = georef((z1=z1, z2=z2))
  w1 = Dict(:z1 => 10, :z2 => 0.1)
  w2 = Dict(:z1 => 0.1, :z2 => 10)
  p1 = partition(𝒮, SLIC(50, 0.001, weights=w1))
  p2 = partition(𝒮, SLIC(50, 0.001, weights=w2))
  @test 50 ≤ length(p1) ≤ 60
  @test 50 ≤ length(p2) ≤ 60
  
  if visualtests
    @test_reference "data/slic-w1.png" plot(p1)
    @test_reference "data/slic-w2.png" plot(p2)
  end

  # test GeoClustering.slic_srecursion function
  k = 20
  l = [10.0, 100.0, 1000.0]
  s = GeoClustering.slic_srecursion(k, l)
  @test s[1] == 10/3 && s[2] == 100/3 && s[3] == 1000/3

  # test SLIC with known problematic domain
  k = 1
  m = 0.000001
  x = LinRange(550350.6224548942, 552307.2106300013, 1200)
  y = LinRange(9.35909841165263e6, 9.36050447440832e6, 1200)
  z = LinRange(-44.90690201082941, 351.4007207008662, 1200)
  Z = (x = x, y = y, z = z, a = rand(1200))
  𝒮 = georef(Z, (:x, :y, :z))
  s = GeoClustering.slic_spacing(𝒮, SLIC(k, m))
  c = GeoClustering.slic_initialization(𝒮, s)
  @test s[1] == 1956.5881751070265 && s[2] == 1406.062755689025 && s[3] == 396.3076227116956
  @test length(c) == 1

  # visual SLIC test for the μCT image
  k = 45
  m = 0.55
  μCT = load("data/muCT.tif")
  𝒮 = georef((μCT = Float64.(μCT),))
  C = cluster(𝒮, SLIC(45, 0.55))

  if visualtests
    @test_reference "data/muCT.png" plot(C)
  end
end
