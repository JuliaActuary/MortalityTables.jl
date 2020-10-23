using MortalityTables
using Test
using JSON
using CSV
using Pkg.Artifacts

# This is the path to the Artifacts.toml we will manipulate
MT_artifact_toml = joinpath(pkgdir(MortalityTables), "Artifacts.toml")
tbl_dir = artifact_hash("SOA_Tables", MT_artifact_toml) |> artifact_path

include("CSV.jl")
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
