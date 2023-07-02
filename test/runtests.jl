using MortalityTables
using Test
using JSON
using CSV
using Pkg.Artifacts

# This is the path to the Artifacts.toml we will manipulate
MT_artifact_toml = joinpath(pkgdir(MortalityTables), "Artifacts.toml")
soa_tbl_dir = artifact_hash("mort.soa.org", MT_artifact_toml) |> artifact_path

include("CSV.jl")
include("basic.jl")
include("XTbML.jl")
include("parameterized_models.jl")
include("distribution.jl")
include("life_expectancy.jl")
include("dukes_macdonald.jl")

# load tables to be used in subsequent tests
@test isa(MortalityTables.table(1), MortalityTables.MortalityTable)

include("projection_scale.jl")
include("MortalityTables.jl")
include("get_SOA_tables.jl")
