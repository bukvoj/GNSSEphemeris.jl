using GNSSEphemeris
using Test

@testset "GNSSEphemeris.jl" begin
    include("lineofsight.jl")
    include("lookangles.jl")
    include("svpos.jl")
end
