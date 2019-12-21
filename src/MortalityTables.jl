module MortalityTables

include("Mortality.jl")


export MortalityTable,
    q, p, qx,
    TableMetaData,
    MortalityAssumption,
    Uniform,Balducci,Constant
end # module
