@testset "face reconstruction" begin
    @test NSCavity.face(Central(), 1.0, 3.0,  5.0) == 2.0
    @test NSCavity.face(Blend(1.0), 1.0, 3.0,  1.0) == 1.0   # upwind side is low side
    @test NSCavity.face(Blend(1.0), 1.0, 3.0, -1.0) == 3.0   # upwind side is high side
    @test NSCavity.face(Blend(0.0), 1.0, 3.0,  1.0) == 2.0   # γ = 0 is central
end

@testset "uniform flow has zero convection" begin
    g = Grid(8)
    U = ones(9, 10); V = zeros(10, 9)
    Cu = zero(U); Cv = zero(V)
    NSCavity.convection!(Cu, Cv, U, V, g, Central())
    @test maximum(abs, Cu) == 0.0
    @test maximum(abs, Cv) == 0.0
end