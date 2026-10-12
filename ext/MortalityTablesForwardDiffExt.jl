module MortalityTablesForwardDiffExt

using MortalityTables, ForwardDiff

# Domain checks use values, including through nested differentiation. The model
# fields and calculations still retain the original dual numbers.
MortalityTables._primal(x::ForwardDiff.Dual) = MortalityTables._primal(ForwardDiff.value(x))

end
