using NSCavity, Plots
const OUT = joinpath(@__DIR__, "..", "results")

function temporal_order(; N = 32, Re = 100.0, T = 0.4, tramp = 0.2,
                          ndt = [100, 200, 400, 800], nref = 3200)
    run(ts, n) = solve_cavity(N = N, Re = Re, tend = T, dt = T / n, tramp = tramp,
                              tscheme = ts).fields.U
    ref = run(AB2(), nref)
    pl = plot(xscale = :log10, yscale = :log10, xlabel = "dt", ylabel = "Linf error of u at T",
              legend = :bottomright)
    for (name, ts) in (("Euler", Euler()), ("AB2", AB2()))
        errs = [maximum(abs, run(ts, n) .- ref) for n in ndt]
        orders = [log(errs[k] / errs[k+1]) / log(2) for k in 1:length(ndt)-1]
        println(name, " errors: ", errs, "\n", name, " observed orders: ", orders)
        plot!(pl, T ./ ndt, errs, marker = :o, label = name)
    end
    dts = T ./ ndt
    plot!(pl, dts, 0.5 .* dts ./ dts[end] .* 1e-2, ls = :dash, label = "order 1")
    plot!(pl, dts, (dts ./ dts[end]).^2 .* 1e-4, ls = :dash, label = "order 2")
    savefig(pl, joinpath(OUT, "temporal_order.png"))
end

temporal_order()