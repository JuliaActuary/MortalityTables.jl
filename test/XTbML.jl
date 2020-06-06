@testset "XTbML" begin
    @testset "dict parse" begin
        g = MortalityTables.get_and_parse
        d = Dict(:a => "1.2", :c => "a")

        @test g(d,:a) == 1.2
        @test ismissing(g(d,:b))
        @test ismissing(g(d,:b))
        @test_throws ArgumentError g(d,:c)

    end
    
    @testset "XTbML loading" begin
        pth = joinpath(
            dirname(pathof(MortalityTables)),
            "tables",
            "SOA",
            "t1076.xml",
        )
        file = MortalityTables.open_and_read(pth) |> MortalityTables.getXML
        xtbl = MortalityTables.parseXTbMLTable(file, pth)
        @test isa(xtbl, MortalityTables.XTbMLTable)
    end

    @testset "Ultimate Only" begin
        pth = joinpath(
        dirname(pathof(MortalityTables)),
        "tables",
        "SOA",
        "t17.xml",
    )
        file = MortalityTables.open_and_read(pth) |> MortalityTables.getXML
        xtbl = MortalityTables.parseXTbMLTable(file, pth)

        mt = MortalityTables.XTbML_Table_To_MortalityTable(xtbl)
        @test isa(mt, MortalityTable)

        @test mt.ultimate[0] ≈ 0.00245
        @test mt.ultimate[100] ≈ 1.0
        @test_throws BoundsError mt.ultimate[101]
    end

    @testset "XTbML to MortalityTable" begin
        @testset "Select and Ultimate" begin
            pth = joinpath(
                dirname(pathof(MortalityTables)),
                "tables",
                "SOA",
                "t1076.xml",
            )
            file = MortalityTables.open_and_read(pth) |> MortalityTables.getXML
            xtbl = MortalityTables.parseXTbMLTable(file, pth)

            mt = MortalityTables.XTbML_Table_To_MortalityTable(xtbl)
            @test isa(mt, MortalityTable)

            @test mt.select[35][35] ≈ 0.00037
            @test mt.ultimate[16] ≈ 0.00041
            @test mt.select[35][59] ≈ 0.00508
            @test mt.select[35][60] ≈ 0.00621 
        end
        
        @testset "Ultimate Only, not begin at age 0" begin
            pth = joinpath(
                dirname(pathof(MortalityTables)),
                "tables",
                "SOA",
                "t18.xml",
            )
            file = MortalityTables.open_and_read(pth) |> MortalityTables.getXML
            xtbl = MortalityTables.parseXTbMLTable(file, pth)
            mt = MortalityTables.XTbML_Table_To_MortalityTable(xtbl)
            @test isa(mt, MortalityTable)
        end
    end
end
