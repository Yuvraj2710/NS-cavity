abstract type TimeScheme end
struct Euler <: TimeScheme end
struct AB2   <: TimeScheme end

# Adams-Bashforth weights for F^n and F^(n-1); the first step falls back to Euler.
ab_coeffs(::Euler, it) = (1.0, 0.0)
ab_coeffs(::AB2,  it) = it == 1 ? (1.0, 0.0) : (1.5, -0.5)

# Conservative dt estimates. The AB2 limits are my estimate, so tune them if it blows up.
stable_dt(g::Grid, p::Params, ::Euler; safety = 0.5) = safety * min(g.h / p.ulid, g.h^2 / (4p.nu), 2p.nu / p.ulid^2)
stable_dt(g::Grid, p::Params, ::AB2; safety = 0.5) = safety * min(g.h / p.ulid, g.h^2 / (8p.nu), p.nu / p.ulid^2)