using Test
using MortalityTables

const mt = MortalityTables

tables = mt.Tables()
cso2001 = tables["2001 CSO Super Preferred Select and Ultimate - Male Nonsmoker, ANB"]
vbt2001 = tables["2001 VBT Residual Standard Select and Ultimate - Male Nonsmoker, ANB"]
cso1980 = tables["1980 CSO Basic Table – Male, ANB"]

@test mt.qx(cso1980,35,1) == .00118
@test mt.qx(cso1980,35,61) == .27302
@test mt.qx(cso1980,95) == .27302
@test ismissing(mt.qx(cso1980,35,95))
@test ismissing(mt.qx(cso1980,101))

@test mt.qx(cso2001,35,1) == .00037
@test mt.qx(cso2001,35,61) == .26719
@test mt.qx(cso2001,16) == .00041
@test mt.qx(cso2001,95) == .26719
@test ismissing(mt.qx(cso2001,15))
@test ismissing(mt.qx(cso2001,150))
@test ismissing(mt.qx(cso2001,35,95))

@test mt.qx(vbt2001,35,1) == .00036
@test mt.qx(vbt2001,35,61) == .24298
@test mt.qx(vbt2001,95) == .24298
@test ismissing(mt.qx(vbt2001,35,95))
@test ismissing(mt.qx(vbt2001,150))

#this is to check trailing whitespace, as the source file has trailing whitespace in it
vbt2001su = tables["2001 VBT Select and Ultimate - Male Nonsmoker, ANB"]
vbt2001su == tables["2001 VBT Select and Ultimate - Male Nonsmoker, ANB"]