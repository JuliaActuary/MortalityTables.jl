module MortalityTables

include("Mortality.jl")


export MortalityTable,
    q, p, qx,
    omega,ω,
    TableMetaData,
    MortalityAssumption,
    Uniform,Balducci,Constant
end # module
