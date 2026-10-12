module MortalityTablesForwardDiffExt

using MortalityTables, ForwardDiff

# Domain checks use values, including through nested differentiation. The model
# fields and calculations still retain the original dual numbers.
MortalityTables._primal(x::ForwardDiff.Dual) = MortalityTables._primal(ForwardDiff.value(x))

# ForwardDiff 0.10's zero comparison ignores partials; a zero mixture weight
# still contributes to a parameter derivative when its seed is nonzero.
MortalityTables._isconstantzero(x::ForwardDiff.Dual) =
    MortalityTables._isconstantzero(ForwardDiff.value(x)) &&
    all(MortalityTables._isconstantzero, ForwardDiff.partials(x))

end
