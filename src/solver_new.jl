function centerlines(U, V, g::Grid)
    N, h = g.N, g.h
    m = N ÷ 2
    return (; yu   = [(j - 0.5) * h for j in 1:N], u_cl = U[m+1, 2:N+1], xv   = [(i - 0.5) * h for i in 1:N], v_cl = V[2:N+1, m+1])
end


function solve_cavity(; N = 64, Re = 100.0, ulid = 1.0, tol = 1e-7, maxit = 400_000, scheme::ConvectionScheme = Central(), tscheme::TimeScheme = AB2(), dt = nothing, tend = nothing, tramp = 0.0, verbose = false)
    g = Grid(N)
    p = Params(Re = Re, ulid = ulid)
    nu = p.nu
    f = Fields(g)
    U, V, P = f.U, f.V, f.P
    ws = PressureWS(g)
    Φ = zero(P)
    Cu = zero(U); Cv = zero(V); Lu = zero(U); Lv = zero(V)
    Fu = zero(U); Fv = zero(V); Fu_prev = zero(U); Fv_prev = zero(V)
    Uold = similar(U)

    dt_ = dt === nothing ? stable_dt(g, p, tscheme) : float(dt)
    transient = tend !== nothing
    nsteps = transient ? round(Int, tend / dt_) : maxit
    lid(t) = tramp > 0 ? ulid * sin(0.5π * min(t / tramp, 1.0))^2 : ulid

    t = 0.0
    apply_bc!(U, V, g, lid(t))
    res = Inf
    it = 0
    while it < nsteps && (transient || res > tol)
        it += 1
        Uold .= U
        convection!(Cu, Cv, U, V, g, scheme)
        laplacian!(Lu, Lv, U, V, g)
        a, b = ab_coeffs(tscheme, it)
        @. Fu = -Cu + nu * Lu
        @. Fv = -Cv + nu * Lv
        U .+= dt_ .* (a .* Fu .+ b .* Fu_prev)
        V .+= dt_ .* (a .* Fv .+ b .* Fv_prev)
        Fu_prev .= Fu; Fv_prev .= Fv
        add_pressure_gradient!(U, V, P, g, -dt_)       # incremental: old pressure gradient
        t += dt_
        apply_bc!(U, V, g, lid(t))
        project!(U, V, Φ, ws, g, dt_, lid(t))
        P .+= Φ
        res = maximum(abs, U .- Uold) / dt_
        verbose && it % 1000 == 0 && println("it=$it  t=$t  residual=$res")
    end

    divergence!(ws.D, U, V, g)
    return (; centerlines(U, V, g)..., iters = it, residual = res, dt = dt_, t = t, divmax = maximum(abs, ws.D), grid = g, fields = f)
end