struct PressureWS
    D::Matrix{Float64}
    b::Matrix{Float64}
    r::Matrix{Float64}
    Ap::Matrix{Float64}
    d::Matrix{Float64}
end

PressureWS(g::Grid) = PressureWS(ntuple(_ -> zeros(g.N, g.N), 5)...)

function divergence!(D, U, V, g::Grid)
    N, h = g.N, g.h
    @inbounds for cj in 1:N, ci in 1:N
        D[ci, cj] = (U[ci+1, cj+1] - U[ci, cj+1] + V[ci+1, cj+1] - V[ci+1, cj]) / h
    end
    return nothing
end

# -Laplacian with homogeneous Neumann BCs (singular: constants are in the null space).
function neglap!(out, P, g::Grid)
    N, h = g.N, g.h
    @inbounds for j in 1:N, i in 1:N
        c = 0.0; s = 0.0
        if i > 1; c += 1; s += P[i-1, j]; end
        if i < N; c += 1; s += P[i+1, j]; end
        if j > 1; c += 1; s += P[i, j-1]; end
        if j < N; c += 1; s += P[i, j+1]; end
        out[i, j] = (c * P[i, j] - s) / h^2
    end
    return nothing
end

# Solves neglap(P) = ws.b by conjugate gradients; returns the iteration count.
function cg!(P, ws::PressureWS, g::Grid; tol = 1e-10, maxit = 5000)
    b, r, Ap, d = ws.b, ws.r, ws.Ap, ws.d
    b .-= sum(b) / length(b)                       # compatibility for Neumann problem
    neglap!(Ap, P, g)
    r .= b .- Ap
    d .= r
    rs = dot(r, r)
    nb = sqrt(dot(b, b)) + eps()
    k = 0
    while k < maxit && sqrt(rs) > tol * nb
        k += 1
        neglap!(Ap, d, g)
        α = rs / dot(d, Ap)
        P .+= α .* d
        r .-= α .* Ap
        rs_new = dot(r, r)
        d .= r .+ (rs_new / rs) .* d
        rs = rs_new
    end
    P .-= sum(P) / length(P)                       # fix the free constant
    return k
end

# Projection: solve lap(p) = div(u*)/dt and make the velocity divergence-free.
function project!(U, V, P, ws::PressureWS, g::Grid, dt, ulid)
    N, h = g.N, g.h
    divergence!(ws.D, U, V, g)
    ws.b .= .-ws.D ./ dt
    iters = cg!(P, ws, g)
    @inbounds for j in 2:N+1, i in 2:N
        U[i, j] -= dt * (P[i, j-1] - P[i-1, j-1]) / h
    end
    @inbounds for j in 2:N, i in 2:N+1
        V[i, j] -= dt * (P[i-1, j] - P[i-1, j-1]) / h
    end
    apply_bc!(U, V, g, ulid)
    return iters
end