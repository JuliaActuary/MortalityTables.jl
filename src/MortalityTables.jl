module MortalityTables
using CSV
using HTTP
using Transducers
using OffsetArrays
using Parsers
using QuadGK
using UnPack
using XMLDict

include("MetaData.jl")
include("death_distribution.jl")
include("MortalityTable.jl")
include("XTbML.jl")
include("CSV.jl")
include("get_SOA_table.jl")
include("parameterized_models.jl")

export MortalityTable,
    survival,
    decrement,
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
    Makeham, Gompertz, MakehamGompertz,
    hazard,cumhazard,
    mortality_vector
end # module
