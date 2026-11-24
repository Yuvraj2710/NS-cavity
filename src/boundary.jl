# Wall BCs: normal velocity is zero on the wall, tangential velocity via ghost cells.
function apply_bc!(U, V, g::Grid, ulid)
    N = g.N
    U[1, :] .= 0.0;  U[N+1, :] .= 0.0
    V[:, 1] .= 0.0;  V[:, N+1] .= 0.0
    @inbounds for i in 1:N+1
        U[i, 1]   = -U[i, 2]                 # bottom wall, u = 0
        U[i, N+2] = 2ulid - U[i, N+1]        # moving lid, u = ulid
    end
    @inbounds for j in 1:N+1
        V[1, j]   = -V[2, j]                 # left wall, v = 0
        V[N+2, j] = -V[N+1, j]               # right wall, v = 0
    end
    return nothing
end