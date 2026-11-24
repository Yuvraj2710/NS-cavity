module NSCavity
export solve_cavity

# Layout (h = 1/N):  U[i,j+1] = u at x=(i-1)h, y=(j-.5)h   (i=1..N+1, j=1..N)
#                    V[i+1,j] = v at x=(i-.5)h, y=(j-1)h   (i=1..N,  j=1..N+1)
#                    P[i,j]   = p at cell centres.  One ghost row/col around U, V.

function apply_bc!(U, V, N, ulid)
    U[1, :] .= 0; U[N+1, :] .= 0; V[:, 1] .= 0; V[:, N+1] .= 0
    for i in 1:N+1
        U[i, 1]   = -U[i, 2]                  # no-slip bottom
        U[i, N+2] = 2ulid - U[i, N+1]         # moving lid
    end
    for j in 1:N+1
        V[1, j]   = -V[2, j]
        V[N+2, j] = -V[N+1, j]
    end
end

function momentum_rhs!(Fu, Fv, U, V, N, h, nu)
    for j in 2:N+1, i in 2:N
        ue = (U[i,j] + U[i+1,j])/2; uw = (U[i-1,j] + U[i,j])/2
        un = (U[i,j] + U[i,j+1])/2; us = (U[i,j-1] + U[i,j])/2
        vn = (V[i,j] + V[i+1,j])/2; vs = (V[i,j-1] + V[i+1,j-1])/2
        conv = (ue^2 - uw^2)/h + (un*vn - us*vs)/h
        diff = (U[i+1,j] + U[i-1,j] + U[i,j+1] + U[i,j-1] - 4U[i,j])/h^2
        Fu[i,j] = -conv + nu*diff
    end
    for j in 2:N, i in 2:N+1
        ue = (U[i,j] + U[i,j+1])/2;   uw = (U[i-1,j] + U[i-1,j+1])/2
        ve = (V[i,j] + V[i+1,j])/2;   vw = (V[i-1,j] + V[i,j])/2
        vn = (V[i,j] + V[i,j+1])/2;   vs = (V[i,j-1] + V[i,j])/2
        conv = (ue*ve - uw*vw)/h + (vn^2 - vs^2)/h
        diff = (V[i+1,j] + V[i-1,j] + V[i,j+1] + V[i,j-1] - 4V[i,j])/h^2
        Fv[i,j] = -conv + nu*diff
    end
end

function divergence!(D, U, V, N, h)
    for cj in 1:N, ci in 1:N
        D[ci,cj] = (U[ci+1,cj+1] - U[ci,cj+1] + V[ci+1,cj+1] - V[ci+1,cj]) / h
    end
end

function neglap!(out, P, N, h)          # -Laplacian, Neumann BC (singular)
    for j in 1:N, i in 1:N
        c = 0.0; s = 0.0
        if i > 1; c += 1; s += P[i-1,j]; end
        if i < N; c += 1; s += P[i+1,j]; end
        if j > 1; c += 1; s += P[i,j-1]; end
        if j < N; c += 1; s += P[i,j+1]; end
        out[i,j] = (c*P[i,j] - s)/h^2
    end
end

function cg!(P, b, N, h; tol=1e-10, maxit=2000)
    b .-= sum(b)/length(b)               # enforce compatibility
    r = similar(b); Ap = similar(b)
    neglap!(Ap, P, N, h); r .= b .- Ap
    d = copy(r); rs = sum(abs2, r); nb = sqrt(sum(abs2, b)) + eps()
    for _ in 1:maxit
        sqrt(rs) < tol*nb && break
        neglap!(Ap, d, N, h)
        a = rs / sum(d .* Ap)
        P .+= a .* d; r .-= a .* Ap
        rs_new = sum(abs2, r)
        d .= r .+ (rs_new/rs) .* d; rs = rs_new
    end
    P .-= sum(P)/length(P)
end

function solve_cavity(; N=64, Re=100.0, ulid=1.0, tol=1e-7, maxit=400_000)
    h = 1/N; nu = ulid/Re
    dt = 0.5*min(h/ulid, h^2/(4nu), 2nu/ulid^2)
    U = zeros(N+1, N+2); V = zeros(N+2, N+1); P = zeros(N, N)
    Fu = zero(U); Fv = zero(V); D = zero(P); Uold = similar(U)
    apply_bc!(U, V, N, ulid)
    res = Inf; it = 0
    while it < maxit && res > tol
        it += 1; Uold .= U
        momentum_rhs!(Fu, Fv, U, V, N, h, nu)
        U .+= dt .* Fu; V .+= dt .* Fv          # predictor (Fu, Fv are zero on boundaries)
        apply_bc!(U, V, N, ulid)
        divergence!(D, U, V, N, h)
        cg!(P, -D ./ dt, N, h)                  # solves  lap(p) = div/dt
        for j in 2:N+1, i in 2:N; U[i,j] -= dt*(P[i,j-1] - P[i-1,j-1])/h; end
        for j in 2:N, i in 2:N+1; V[i,j] -= dt*(P[i-1,j] - P[i-1,j-1])/h; end
        apply_bc!(U, V, N, ulid)
        res = maximum(abs.(U .- Uold)) / dt
    end
    m = N ÷ 2                                   # N must be even
    return (; yu = [(j-0.5)*h for j in 1:N], u_cl = U[m+1, 2:N+1],
              xv = [(i-0.5)*h for i in 1:N], v_cl = V[2:N+1, m+1],
              iters = it, residual = res, dt = dt)
end
end