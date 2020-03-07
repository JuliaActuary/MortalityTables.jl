@testset "XTbML" begin

    @testset "XTbML loading" begin
        pth = joinpath(dirname(pathof(MortalityTables)), "tables", "SOA","t1076.xml")
        file = MortalityTables.open_and_read(pth) |> MortalityTables.getXML
        xtbl = MortalityTables.parseXTbMLTable(file,pth)
        @test isa(xtbl,MortalityTables.XTbMLTable)
    end


    @testset "XTbML rates" begin
        pth = joinpath(dirname(pathof(MortalityTables)), "tables", "SOA","t1076.xml")
        file = MortalityTables.open_and_read(pth) |> MortalityTables.getXML
        xtbl = MortalityTables.parseXTbMLTable(file,pth)

        @test MortalityTables.q_select(xtbl,35,1) ≈ .00037
        @test MortalityTables.q_ultimate(xtbl,16) ≈ .00041
        @test MortalityTables.q_select(xtbl,35,25) ≈ 0.00508
        @test MortalityTables.q_select(xtbl,35,26) ≈ 	0.00621
    end

    @testset "XTbML to MortalityTable" begin
        @testset "Select and Ultimate" begin
            pth = joinpath(dirname(pathof(MortalityTables)), "tables", "SOA","t1076.xml")
            file = MortalityTables.open_and_read(pth) |> MortalityTables.getXML
            xtbl = MortalityTables.parseXTbMLTable(file,pth)

            mt = MortalityTables.XTbML_Table_To_MortalityTable(xtbl)
            @test isa(mt,MortalityTable)
        end
        @testset "Ultimate Only" begin
            pth = joinpath(dirname(pathof(MortalityTables)), "tables", "SOA","t17.xml")
            file = MortalityTables.open_and_read(pth) |> MortalityTables.getXML
            xtbl = MortalityTables.parseXTbMLTable(file,pth)

            mt = MortalityTables.XTbML_Table_To_MortalityTable(xtbl)
            @test isa(mt,MortalityTable)
        end
        @testset "Ultimate Only, not begin at age 0" begin
            pth = joinpath(dirname(pathof(MortalityTables)), "tables", "SOA","t18.xml")
            file = MortalityTables.open_and_read(pth) |> MortalityTables.getXML
            xtbl = MortalityTables.parseXTbMLTable(file,pth)
            mt = MortalityTables.XTbML_Table_To_MortalityTable(xtbl)
            @test isa(mt,MortalityTable)
        end
    end
end
