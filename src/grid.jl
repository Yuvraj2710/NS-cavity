# Staggered (MAC) layout on the unit square, h = 1/N, one ghost layer on U and V:
#   U[i, j+1] : u at x=(i-1)h,   y=(j-1/2)h   i=1..N+1, j=1..N   (U is (N+1)x(N+2))
#   V[i+1, j] : v at x=(i-1/2)h, y=(j-1)h     i=1..N,   j=1..N+1 (V is (N+2)x(N+1))
#   P[i, j]   : p at cell centres                                 (P is NxN)

struct Grid
    N::Int
    h::Float64
end

function Grid(N::Int)
    iseven(N) || throw(ArgumentError("N must be even so the centerline lies on grid lines"))
    return Grid(N, 1.0 / N)
end

struct Params
    Re::Float64
    ulid::Float64
    nu::Float64
end

Params(; Re::Real = 100.0, ulid::Real = 1.0) = Params(Re, ulid, ulid / Re)  # L = 1

struct Fields
    U::Matrix{Float64}
    V::Matrix{Float64}
    P::Matrix{Float64}
end

Fields(g::Grid) = Fields(zeros(g.N + 1, g.N + 2), zeros(g.N + 2, g.N + 1), zeros(g.N, g.N))