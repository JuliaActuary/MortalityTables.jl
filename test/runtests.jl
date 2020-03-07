using MortalityTables
using Test

include("basic.jl")
include("XTbML.jl")

#load tables to be used in subsequent tests
tables = MortalityTables.tables()
@test length(tables) > 0

include("MortalityTables.jl")
