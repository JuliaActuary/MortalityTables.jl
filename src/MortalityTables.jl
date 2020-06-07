module MortalityTables
using HTTP
using Transducers
using OffsetArrays
using QuadGK
using XMLDict

include("MetaData.jl")
include("death_distribution.jl")
include("MortalityTable.jl")
include("XTbML.jl")
include("get_SOA_table.jl")
include("parameterized_models.jl")

export MortalityTable,
    survivorship,
    cumulative_decrement,
    omega,
    TableMetaData,
    SelectMortality,
    UltimateMortality,
    Balducci,
    Uniform,
    Constant,
    DeathDistribution,
    get_SOA_table,
    get_SOA_table!,
    Makeham, Gompertz
end # module
