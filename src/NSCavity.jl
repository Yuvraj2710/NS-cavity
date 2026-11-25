module NSCavity
export ConvectionScheme, Central, Blend
export TimeScheme, Euler, AB2
using LinearAlgebra

include("grid.jl")
include("boundary.jl")
include("convection.jl")
include("diffusion.jl")
include("pressure.jl")
include("timestepping.jl")
include("solver.jl")

export Grid, Params, Fields
export ConvectionScheme, Central, Blend
export solve_cavity

end