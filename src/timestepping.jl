# Conservative step for explicit Euler: convective CFL, diffusive limit, and the
# 2*nu/u^2 stability bound of forward Euler with central convection.
stable_dt(g::Grid, p::Params; safety = 0.5) =
    safety * min(g.h / p.ulid, g.h^2 / (4p.nu), 2p.nu / p.ulid^2)