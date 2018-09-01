import Test
import MortalityTables

cso_path = joinpath(raw"C:\Users\alecl\AppData\Local\Julia-1.0.0\MortalityTables\tables\t1076.xml")
cso = MortalityTables.loadXTbMLTable(cso_path)



Test.@test MortalityTables.qx(cso,35,1) == .00037
Test.@test MortalityTables.qx(cso,35,61) == .26719
Test.@test MortalityTables.qx(cso,95) == .26719
Test.@test ismissing(MortalityTables.qx(cso,35,95))
