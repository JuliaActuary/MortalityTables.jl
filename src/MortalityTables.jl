module MortalityTables
using HTTP
using Transducers
using OffsetArrays
using Parsers
using QuadGK
using Requires
using UnPack
using XMLDict

include("MetaData.jl")
include("death_distribution.jl")
include("MortalityTable.jl")
include("XTbML.jl")
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

# lazy load part of the package
function __init__()
    @require CSV="336ed68f-0bac-5ca0-87d4-7b16caf5d00b" include("CSV.jl")
end

end # module
