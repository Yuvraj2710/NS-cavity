using NSCavity, Plots
include("cavity_ghia.jl")   # reuses interp1 and profiles (also reruns its Re=100 comparison)

# Self-convergence: Ghia's table is only accurate to ~1e-3, so it cannot show
# convergence order. Compare each N against a finer reference at the same points.
function self_convergence(Re; Ns = [16, 32, 64], Nref = 128)
    ref = profiles(solve_cavity(N = Nref, Re = Float64(Re)))
    pts = range(0.05, 0.95, length = 19)
    errs = map(Ns) do N
        p = profiles(solve_cavity(N = N, Re = Float64(Re)))
        maximum(abs(interp1(p.ys, p.us, y) - interp1(ref.ys, ref.us, y)) for y in pts)
    end
    hs = 1.0 ./ Ns
    orders = [log(errs[k] / errs[k+1]) / log(2) for k in 1:length(Ns)-1]
    println("N, error: ", collect(zip(Ns, errs)), "\nobserved orders: ", orders)
    pl = plot(hs, errs, xscale = :log10, yscale = :log10, marker = :o,
              xlabel = "h", ylabel = "Linf error of u vs N=$Nref", label = "solver")
    plot!(pl, hs, errs[1] .* (hs ./ hs[1]) .^ 2, ls = :dash, label = "order 2")
    savefig(pl, joinpath(OUT, "convergence_re$(Re).png"))
end

self_convergence(100)