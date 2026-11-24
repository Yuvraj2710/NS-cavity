# Five-point Laplacians of u and v (viscosity is applied by the caller).
function laplacian!(Lu, Lv, U, V, g::Grid)
    N, h = g.N, g.h
    ih2 = 1.0 / h^2
    @inbounds for j in 2:N+1, i in 2:N
        Lu[i, j] = (U[i+1, j] + U[i-1, j] + U[i, j+1] + U[i, j-1] - 4U[i, j]) * ih2
    end
    @inbounds for j in 2:N, i in 2:N+1
        Lv[i, j] = (V[i+1, j] + V[i-1, j] + V[i, j+1] + V[i, j-1] - 4V[i, j]) * ih2
    end
    return nothing
end