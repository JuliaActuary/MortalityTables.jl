module MortalityTables
using Memoize
using OffsetArrays
using Parsers
using QuadGK
using Requires
import StringDistances
using UnPack
using XMLDict
using Pkg.Artifacts

include("table_source_map.jl")
include("MetaData.jl")
include("death_distribution.jl")
include("MortalityTable.jl")
include("dukes_macdonald.jl")
include("XTbML.jl")
include("get_SOA_table.jl")
include("parameterized_models.jl")
include("life_expectancy.jl")

table_dirs = Dict(
    "mort.soa.org" => artifact"mort.soa.org",
)

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
    Makeham, Gompertz, MakehamGompertz,
    hazard, cumhazard,
    mortality_vector

# lazy load part of the package
function __init__()
    @require CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b" include("CSV.jl")
end

end # module
