import Test
import MortalityTables

tables = MortalityTables.Tables()
cso = tables["2001 CSO Super Preferred Select and Ultimate - Male Nonsmoker, ANB"]


Test.@test MortalityTables.qx(cso,35,1) == .00037
Test.@test MortalityTables.qx(cso,35,61) == .26719
Test.@test MortalityTables.qx(cso,95) == .26719
Test.@test ismissing(MortalityTables.qx(cso,35,95))
