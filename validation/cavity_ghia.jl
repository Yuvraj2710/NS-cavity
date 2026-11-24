using NSCavity, CSV, DataFrames, Plots

const DATA = joinpath(@__DIR__, "data")
const OUT  = joinpath(@__DIR__, "..", "results")

function interp1(xs, ys, x)                      # linear interpolation, xs ascending
    k = clamp(searchsortedlast(xs, x), 1, length(xs) - 1)
    t = (x - xs[k]) / (xs[k+1] - xs[k])
    return (1 - t) * ys[k] + t * ys[k+1]
end

# Add wall values so Ghia points at y,x = 0 and 1 are inside the interpolation range.
profiles(s) = (ys = [0.0; s.yu; 1.0], us = [0.0; s.u_cl; 1.0],
               xs = [0.0; s.xv; 1.0], vs = [0.0; s.v_cl; 0.0])

function compare(Re::Int, Ns; scheme = Central())
    gu = CSV.read(joinpath(DATA, "ghia_re$(Re)_u.csv"), DataFrame)   # columns: y,u
    gv = CSV.read(joinpath(DATA, "ghia_re$(Re)_v.csv"), DataFrame)   # columns: x,v
    sols = Dict{Int,Any}()
    for N in Ns
        @info "Solving" Re N
        sols[N] = solve_cavity(N = N, Re = Float64(Re), scheme = scheme)
    end
    println("N     Linf(u)      Linf(v)")
    for N in Ns
        p = profiles(sols[N])
        eu = maximum(abs(interp1(p.ys, p.us, y) - u) for (y, u) in zip(gu.y, gu.u))
        ev = maximum(abs(interp1(p.xs, p.vs, x) - v) for (x, v) in zip(gv.x, gv.v))
        println(rpad(N, 6), rpad(round(eu, sigdigits = 3), 13), round(ev, sigdigits = 3))
    end
    Nf = maximum(Ns); p = profiles(sols[Nf])
    pl1 = plot(p.us, p.ys, label = "N=$Nf", xlabel = "u", ylabel = "y", title = "u at x=0.5")
    scatter!(pl1, gu.u, gu.y, label = "Ghia et al. (1982)")
    pl2 = plot(p.xs, p.vs, label = "N=$Nf", xlabel = "x", ylabel = "v", title = "v at y=0.5")
    scatter!(pl2, gv.x, gv.v, label = "Ghia et al. (1982)")
    savefig(plot(pl1, pl2, layout = (1, 2), size = (900, 400)),
            joinpath(OUT, "ghia_re$(Re).png"))
end

compare(100, [16, 32, 64, 128])