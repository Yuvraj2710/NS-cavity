function centerlines(U, V, g::Grid)
    N, h = g.N, g.h
    m = N ÷ 2
    return (; yu   = [(j - 0.5) * h for j in 1:N],
              u_cl = U[m+1, 2:N+1],              # u along x = 0.5
              xv   = [(i - 0.5) * h for i in 1:N],
              v_cl = V[2:N+1, m+1])              # v along y = 0.5
end

function solve_cavity(; N = 64, Re = 100.0, ulid = 1.0, tol = 1e-7, maxit = 400_000,
                        scheme::ConvectionScheme = Central(), verbose = false)
    g = Grid(N)
    p = Params(Re = Re, ulid = ulid)
    nu = p.nu
    f = Fields(g)
    U, V, P = f.U, f.V, f.P
    ws = PressureWS(g)
    Cu = zero(U); Cv = zero(V); Lu = zero(U); Lv = zero(V)
    Fu = zero(U); Fv = zero(V); Uold = similar(U)

    dt = stable_dt(g, p)
    apply_bc!(U, V, g, p.ulid)
    res = Inf
    it = 0
    while it < maxit && res > tol
        it += 1
        Uold .= U
        convection!(Cu, Cv, U, V, g, scheme)
        laplacian!(Lu, Lv, U, V, g)
        @. Fu = -Cu + nu * Lu
        @. Fv = -Cv + nu * Lv
        U .+= dt .* Fu
        V .+= dt .* Fv
        apply_bc!(U, V, g, p.ulid)
        project!(U, V, P, ws, g, dt, p.ulid)
        res = maximum(abs, U .- Uold) / dt
        verbose && it % 1000 == 0 && println("it=$it  residual=$res")
    end

    divergence!(ws.D, U, V, g)
    return (; centerlines(U, V, g)..., iters = it, residual = res, dt = dt,
              divmax = maximum(abs, ws.D), grid = g, fields = f)
end