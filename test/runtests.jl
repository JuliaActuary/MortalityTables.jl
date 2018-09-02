import Test
import MortalityTables

tables = MortalityTables.Tables()
cso = tables["2001 CSO Super Preferred Select and Ultimate - Male Nonsmoker, ANB"]
vbt = tables["2001 VBT Residual Standard Select and Ultimate - Male Nonsmoker, ANB"]

Test.@test MortalityTables.qx(cso,35,1) == .00037
Test.@test MortalityTables.qx(cso,35,61) == .26719
Test.@test MortalityTables.qx(cso,95) == .26719
Test.@test ismissing(MortalityTables.qx(cso,35,95))

Test.@test MortalityTables.qx(vbt,35,1) == .00036
Test.@test MortalityTables.qx(vbt,35,61) == .24298
Test.@test MortalityTables.qx(vbt,95) == .24298
Test.@test ismissing(MortalityTables.qx(vbt,35,95))
