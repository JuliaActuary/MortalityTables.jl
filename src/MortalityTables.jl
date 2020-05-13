module MortalityTables
using EzXML
include("Mortality.jl")
include("death_distribution.jl")
include("get_SOA_table.jl")

export MortalityTable,
    q,
    p,
    omega,
    ω,
    TableMetaData,
    SelectMortality,
    UltimateMortality,
    MortalityVector,
    MortalityTable,
    Balducci,
    Uniform,
    Constant,
    DeathDistribution,
    get_SOA_table,
    get_SOA_table!,
    eztbl
end # module
