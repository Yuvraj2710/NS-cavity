@testset "lid-driven cavity, Re = 100, N = 32" begin
    s = solve_cavity(N = 32, Re = 100.0, tol = 1e-6)
    @test s.residual <= 1e-6                  # reached steady state
    @test s.divmax < 1e-6                     # projection works
    @test 0.7 < s.u_cl[end] < 1.0             # near-lid velocity
    @test -0.3 < minimum(s.u_cl) < -0.15      # return flow; loose bounds, tighten after Ghia check
end