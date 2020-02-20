using MortalityTables
using Test

@testset "basic `Ultimate` MortalityTable" begin
    tbl = MortalityTable([i/100 for i in 1:100],0)
    @test q(tbl,0,1) ≈ 0.01
    @test q(tbl,0,2) ≈ 1.0 - (1.0 - 0.01) * (1.0 - 0.02)
    @test q(tbl,0) ≈ 0.01

    @test p(tbl,0) ≈ 0.99
    @test p(tbl,0,1) ≈ 0.99
    @test p(tbl,0,2) ≈ 1.0 - q(tbl,0,2)
    @test q(tbl,0,1,1) ≈ 0.02
    @test p(tbl,0,1,1) ≈ 0.98
    @test_throws BoundsError q(tbl,200)
end

@testset "basic `Ultimate` MortalityTable with non-zero start age" begin
    tbl = MortalityTable([i/100 for i in 5:100],5)
    @test ismissing(q(tbl,0))
    @test q(tbl,5) ≈ 0.05
    @test q(tbl,5,1) ≈ 0.05
    @test q(tbl,5,2) ≈ 1.0 - (1.0 - 0.05) * (1.0 - 0.06)


    @test p(tbl,5) ≈ 0.95
    @test p(tbl,5,1) ≈ 0.95
    @test p(tbl,5,2) ≈ 1.0 - q(tbl,5,2)

    @test q(tbl,5,1,1) ≈ 0.06
    @test p(tbl,5,1,1) ≈ 0.94
end

@testset "test table names" begin
    tbl = MortalityTable([i/100 for i in 1:100],0,TableMetaData(name="test"))
    @test q(tbl,0,1) ≈ 0.01
    @test tbl.d.name == "test"
end




@testset "XtbML constructors" begin
    pth = joinpath(dirname(pathof(MortalityTables)), "tables", "SOA","t1076.xml")
    file = MortalityTables.open_and_read(pth) |> MortalityTables.getXML
    xtbl = MortalityTables.parseXTbMLTable(file,pth)
    @test qx(xtbl,35,1) ≈ .00037
    @test qx(xtbl,16) ≈ .00041
    @test qx(xtbl,35,25) ≈ 0.00508
    @test qx(xtbl,35,26) ≈ 	0.00621
end

tables = MortalityTables.tables()

@testset "TableMetaData" begin
    @testset "whitespace managment" begin
        # "2017 Loaded CSO Smoker Distinct Nonsmoker Female ALB" comes with a trailing whitespace
        @test "2017 Loaded CSO Smoker Distinct Nonsmoker Female ALB" in keys(tables)
    end
end

@testset "SOA tables" begin


    cso2001 = tables["2001 CSO Super Preferred Select and Ultimate - Male Nonsmoker, ANB"]
    vbt2001 = tables["2001 VBT Residual Standard Select and Ultimate - Male Nonsmoker, ANB"]
    cso1980 = tables["1980 CSO Basic Table – Male, ANB"]

    @test q(cso1980,35,0,1) ≈ .00118
    @test qx(cso1980,35,1) ≈ .00118
    @test qx(cso1980,35,61) ≈ .27302
    @test qx(cso1980,95) ≈ .27302
    @test ismissing(qx(cso1980,35,95))
    @test ismissing(qx(cso1980,101))
    @test omega(cso1980) == 100
    @test ω(cso1980) == 100


    @test qx(cso2001,35,1) ≈ .00037
    @test qx(cso2001,35,61) ≈ .26719
    @test qx(cso2001,16) ≈ .00041
    @test qx(cso2001,95) ≈ .26719
    @test qx(cso2001,120) ≈ 1.0
    @test ismissing(qx(cso2001,15))
    @test_throws BoundsError qx(cso2001,150)
    @test ismissing(qx(cso2001,35,95))
    @test omega(cso2001) == 120
    @test ω(cso2001) == 120

    @test qx(vbt2001,35,1) ≈ .00036
    @test qx(vbt2001,35,61) ≈ .24298
    @test qx(vbt2001,95) ≈ .24298
    @test qx(vbt2001,120) ≈ 1.0
    @test ismissing(qx(vbt2001,35,95))
    @test_throws BoundsError qx(vbt2001,150)
    @test omega(vbt2001) == 120
    @test ω(vbt2001) == 120

    @test cso2001[29,1:3] == [0.00029, 0.00035, 0.0004]

    @info "qx with range"
    #select only
    @test qx(cso2001,35,1:25) == [0.00037,0.00043,0.00049,0.00057,0.00063,0.0007,0.00077,0.00084,0.00092,0.00101,0.00114,0.00127,0.00143,0.00159,0.00174,0.00188,0.00208,0.00231,0.00251,0.00279,0.00315,0.00357,0.00406,0.00461,0.00508]
    # crosses into ultimate
    # crosses into ultimate
    @test qx(cso2001,35,21:30) == [0.00315,0.00357,0.00406,0.00461,0.00508,0.00621,0.0069,0.00773,0.00867,0.00965]
    @test qx(cso2001,35,1:30) == [0.00037,0.00043,0.00049,0.00057,0.00063,0.0007,0.00077,0.00084,0.00092,0.00101,0.00114,0.00127,0.00143,0.00159,0.00174,0.00188,0.00208,0.00231,0.00251,0.00279,0.00315,0.00357,0.00406,0.00461,0.00508,0.00621,0.0069,0.00773,0.00867,0.00965]


end

@testset "Mortality Assumption" begin

    vbt2001 = tables["2001 VBT Residual Standard Select and Ultimate - Male Nonsmoker, ANB"]


    # q(b,20,0.5) ≈ 0.00077 * # Balducci

    # uniform
    u = MortalityAssumption(vbt2001,Uniform())
    @test q(u,25,1) ≈ 0.00101
    @test p(u,25,0.5) ≈ 1 - 0.00101 * 0.5
    @test p(u,25,1.5) ≈ (1 - 0.00101) * (1 - 0.00104 * 0.5)

    #constant
    c = MortalityAssumption(vbt2001,Constant())
    @test q(c,25,1) ≈ 0.00101
    @test p(c,25,0.5) ≈ (1 - 0.00101) ^ 0.5
    @test p(c,25,1.5) ≈ (1 - 0.00101) * (1 - 0.00104) ^ 0.5

end
