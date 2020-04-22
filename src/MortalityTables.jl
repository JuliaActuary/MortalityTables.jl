module MortalityTables

include("Mortality.jl")
include("death_distribution.jl")


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
    DeathDistribution
end # module
