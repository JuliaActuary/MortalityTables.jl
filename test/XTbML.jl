@testset "XTbML" begin
    @testset "dict parse" begin
        g = MortalityTables.get_and_parse
        d = Dict(:a => "1.2", :c => "a")

        @test g(d, :a) == 1.2
        @test ismissing(g(d, :b))
        @test ismissing(g(d, :b))
        @test_throws Exception g(d, :c)

    end
    
    @testset "XTbML loading" begin
        pth = joinpath(soa_tbl_dir,"t1076.xml")
        file = MortalityTables.open_and_read(pth) |> MortalityTables.getXML
        xtbl = MortalityTables.parseXTbMLTable(file, pth)
        @test isa(xtbl, MortalityTables.XTbMLTable)
    end

    @testset "readXTbML cache" begin
        pth = joinpath(soa_tbl_dir,"t1076.xml")
        # repeated reads of the same path return the identical object
        @test MortalityTables.readXTbML(pth) === MortalityTables.readXTbML(pth)
    end

    @testset "Ultimate Only" begin
        pth = joinpath(soa_tbl_dir,"t17.xml")
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
            pth = joinpath(soa_tbl_dir,"t1076.xml")
            file = MortalityTables.open_and_read(pth) |> MortalityTables.getXML
            xtbl = MortalityTables.parseXTbMLTable(file, pth)

            mt = MortalityTables.XTbML_Table_To_MortalityTable(xtbl)
            @test isa(mt, MortalityTable)

            @test mt.select[35][35] ≈ 0.00037
            @test mt.ultimate[16] ≈ 0.00041
            @test mt.select[35][59] ≈ 0.00508
            @test mt.select[35][60] ≈ 0.00621 
        end
        
        @testset "rates are placed by their labels" begin
            md = MortalityTables.TableMetaData(name = "probe")
            ult = [(age = a, rate = 0.2) for a in 40:50]
            row(durs) = [(issue_age = 40, rates = [(duration = d, rate = d / 100) for d in durs])]
            build(sel, u = ult) = MortalityTables.XTbML_Table_To_MortalityTable(MortalityTables.XTbMLTable(sel, u, md))
            # consecutive durations run into the ultimate rates at the next attained age
            mt = build(row(1:3))
            @test mt.select[40][40:44] == [0.01, 0.02, 0.03, 0.2, 0.2]
            @test eltype(mt.select[40]) == Float64
            # a missing leading or interior duration keeps every other rate at its own age
            @test isequal(build(row(2:3)).select[40][40:43], [missing, 0.02, 0.03, 0.2])
            @test isequal(build(row([1, 3])).select[40][40:43], [0.01, missing, 0.03, 0.2])
            # labels, not positions, place the rates
            @test build(row([2, 1])).select[40][40:42] == [0.01, 0.02, 0.2]
            # a repeated label is ambiguous
            @test_throws "repeat" build(row([1, 2, 2]))
            # grouped issue ages: each row at its own issue age, none in between
            grouped = [(issue_age = a, rates = [(duration = 1, rate = a / 1000)]) for a in (40, 42)]
            mt = build(grouped)
            @test mt.select[40][40] == 0.04 && mt.select[42][42] == 0.042
            @test ismissing(mt.select[41])
            # grouped ultimate ages
            mt = build(nothing, [(age = a, rate = a / 1000) for a in (40, 42)])
            @test mt.ultimate[40] == 0.04 && mt.ultimate[42] == 0.042
            @test ismissing(mt.ultimate[41])
        end

        @testset "Ultimate Only, not begin at age 0" begin
            pth = joinpath(soa_tbl_dir,"t18.xml")
            file = MortalityTables.open_and_read(pth) |> MortalityTables.getXML
            xtbl = MortalityTables.parseXTbMLTable(file, pth)
            mt = MortalityTables.XTbML_Table_To_MortalityTable(xtbl)
            @test isa(mt, MortalityTable)
        end
    end
end
