@testset "Adams-Bashforth weights" begin
    @test NSCavity.ab_coeffs(Euler(), 5) == (1.0, 0.0)
    @test NSCavity.ab_coeffs(AB2(), 1)   == (1.0, 0.0)
    @test NSCavity.ab_coeffs(AB2(), 2)   == (1.5, -0.5)
end

@testset "steady state does not depend on the time scheme" begin
    a = solve_cavity(N = 24, Re = 100.0, tol = 1e-7, tscheme = Euler())
    b = solve_cavity(N = 24, Re = 100.0, tol = 1e-7, tscheme = AB2())
    @test maximum(abs, a.u_cl .- b.u_cl) < 2e-3      # loose; tighten once you see the real gap
    @test a.divmax < 1e-6 && b.divmax < 1e-6
end