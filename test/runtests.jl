import Test
import MortalityTables

tables = MortalityTables.Tables()
cso2001 = tables["2001 CSO Super Preferred Select and Ultimate - Male Nonsmoker, ANB"]
vbt2001 = tables["2001 VBT Residual Standard Select and Ultimate - Male Nonsmoker, ANB"]
cso1980 = tables["1980 CSO Basic Table – Male, ANB"]

Test.@test MortalityTables.qx(cso1980,35,1) == .00118
Test.@test MortalityTables.qx(cso1980,35,61) == .27302
Test.@test MortalityTables.qx(cso1980,95) == .27302
Test.@test ismissing(MortalityTables.qx(cso1980,35,95))
Test.@test ismissing(MortalityTables.qx(cso1980,101))

Test.@test MortalityTables.qx(cso2001,35,1) == .00037
Test.@test MortalityTables.qx(cso2001,35,61) == .26719
Test.@test MortalityTables.qx(cso2001,95) == .26719
Test.@test ismissing(MortalityTables.qx(cso2001,35,95))

Test.@test MortalityTables.qx(vbt2001,35,1) == .00036
Test.@test MortalityTables.qx(vbt2001,35,61) == .24298
Test.@test MortalityTables.qx(vbt2001,95) == .24298
Test.@test ismissing(MortalityTables.qx(vbt2001,35,95))
