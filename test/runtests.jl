using MortalityTables
using Test
using JSON

include("CSVs.jl")
include("basic.jl")
include("XTbML.jl")
include("parameterized_models.jl")
include("distribution.jl")

# load tables to be used in subsequent tests
tables = MortalityTables.tables()
@test length(tables) > 0

include("projection_scale.jl")
include("MortalityTables.jl")
include("get_SOA_tables.jl")
