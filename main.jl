using CSV
using DataFrames
using Random
using StatsBase
using Distributions
using Profile

include("./src/types.jl")
include("./src/init.jl")
include("./src/run_generations.jl")
include("./src/dispersal.jl")
include("./src/selection.jl")
include("./src/reproduction.jl")
include("./src/logging.jl")

function main()
    mpfm::MPFMInstance = init()
    #printsum(mpfm)
    @profile run_mpfm(mpfm)
    Profile.print()
end

main()
