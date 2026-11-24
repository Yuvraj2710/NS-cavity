using Test
using NSCavity

@testset "NSCavity" begin
    include("test_pressure.jl")
    include("test_convection.jl")
    include("test_cavity.jl")
end