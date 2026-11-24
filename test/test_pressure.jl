using Random

@testset "CG recovers a known pressure field" begin
    g = NSCavity.Grid(16)
    ws = NSCavity.PressureWS(g)
    rng = MersenneTwister(1)
    P0 = randn(rng, 16, 16)
    P0 .-= sum(P0) / length(P0)
    NSCavity.neglap!(ws.b, P0, g)
    P = zeros(16, 16)
    NSCavity.cg!(P, ws, g; tol = 1e-12, maxit = 10_000)
    @test maximum(abs, P .- P0) < 1e-8
end