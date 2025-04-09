module GNSSEphemeris
    using LinearAlgebra
    using Dates, TimesDates
    using Geodesy
    using RinexRead

    include("utils.jl")
    include("lineofsight.jl")
    include("lookangles.jl")
    include("ephemeris.jl")

    export getsvpos, lineofsight, lookangles
end
