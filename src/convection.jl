abstract type ConvectionScheme end

struct Central <: ConvectionScheme end

# γ = 0: central, γ = 1: first-order upwind (donor cell). Intermediate: Hirt-Cook blend.
struct Blend <: ConvectionScheme
    γ::Float64
    function Blend(γ)
        0.0 <= γ <= 1.0 || throw(ArgumentError("γ must be in [0, 1]"))
        return new(γ)
    end
end

# Value of the transported quantity at a face between cells with values a (low side)
# and b (high side); vel is the transporting velocity at that face.
@inline face(::Central, a, b, vel) = 0.5 * (a + b)
@inline face(s::Blend, a, b, vel) = 0.5 * (a + b) - 0.5 * s.γ * sign(vel) * (b - a)

# Conservative convective terms  d(uu)/dx + d(vu)/dy  and  d(uv)/dx + d(vv)/dy.
function convection!(Cu, Cv, U, V, g::Grid, s::ConvectionScheme)
    N, h = g.N, g.h
    @inbounds for j in 2:N+1, i in 2:N
        ue = 0.5 * (U[i, j] + U[i+1, j]);   uw = 0.5 * (U[i-1, j] + U[i, j])
        vn = 0.5 * (V[i, j] + V[i+1, j]);   vs = 0.5 * (V[i, j-1] + V[i+1, j-1])
        fe = ue * face(s, U[i, j],   U[i+1, j], ue)
        fw = uw * face(s, U[i-1, j], U[i, j],   uw)
        fn = vn * face(s, U[i, j],   U[i, j+1], vn)
        fs = vs * face(s, U[i, j-1], U[i, j],   vs)
        Cu[i, j] = (fe - fw + fn - fs) / h
    end
    @inbounds for j in 2:N, i in 2:N+1
        ue = 0.5 * (U[i, j] + U[i, j+1]);   uw = 0.5 * (U[i-1, j] + U[i-1, j+1])
        vn = 0.5 * (V[i, j] + V[i, j+1]);   vs = 0.5 * (V[i, j-1] + V[i, j])
        fe = ue * face(s, V[i, j],   V[i+1, j], ue)
        fw = uw * face(s, V[i-1, j], V[i, j],   uw)
        fn = vn * face(s, V[i, j],   V[i, j+1], vn)
        fs = vs * face(s, V[i, j-1], V[i, j],   vs)
        Cv[i, j] = (fe - fw + fn - fs) / h
    end
    return nothing
end