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
include("age_basis.jl")
include("dukes_macdonald.jl")
include("labeled_tables.jl")
include("XTbML.jl")
include("get_SOA_table.jl")
include("parametric_interface.jl")
include("parametric_kernels.jl")
include("parameterized_models.jl")
include("ratio_laws.jl")
include("life_expectancy.jl")

export MortalityTable,
    survival,
    decrement,
    curtate_life_expectancy,
    complete_life_expectancy,
    omega,
    TableMetaData,
    SelectMortality,
    UltimateMortality,
    Balducci,
    UniformDeaths,
    ConstantForce,
    DeathDistribution,
    get_SOA_table,
    Makeham, Gompertz,
    hazard, cumhazard,
    mortality_vector,
    age_nearest_to_age_last,
    age_last_to_age_nearest

end # module
