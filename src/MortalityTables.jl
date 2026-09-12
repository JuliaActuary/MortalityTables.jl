module MortalityTables
using OffsetArrays
using QuadGK
import StringDistances
using XML
using Artifacts

include("table_source_map.jl")
include("MetaData.jl")
include("death_distribution.jl")
include("MortalityTable.jl")
include("dukes_macdonald.jl")
include("XTbML.jl")
include("get_SOA_table.jl")
include("parametric_interface.jl")
include("parameterized_models.jl")
include("life_expectancy.jl")

export MortalityTable,
    survival,
    decrement,
    life_expectancy,
    omega,
    TableMetaData,
    SelectMortality,
    UltimateMortality,
    Balducci,
    Uniform,
    Constant,
    DeathDistribution,
    get_SOA_table,
    Makeham, Gompertz,
    hazard, cumhazard,
    mortality_vector

end # module
