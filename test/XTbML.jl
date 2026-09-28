using XML: XML

@testset "XTbML" begin
    @testset "node helpers" begin
        node(s) = only(XML.elements(XML.parse(s, XML.Node)))

        # a rate cell is a Float64, an empty cell is missing, anything else is an error
        @test MortalityTables._rate(node("<Y t=\"1\">0.5</Y>")) == 0.5
        @test MortalityTables._rate(node("<Y t=\"1\"></Y>")) === missing
        @test MortalityTables._rate(node("<Y t=\"1\"> </Y>")) === missing
        @test_throws ArgumentError MortalityTables._rate(node("<Y t=\"1\">abc</Y>"))

        # a metadata field is nothing when absent, "" when empty, and may carry attributes
        parent = node("<md><a tc=\"4\">x</a><b></b></md>")
        @test MortalityTables._text(parent, "zzz") === nothing
        @test MortalityTables._text(parent, "b") == ""
        @test MortalityTables._text(parent, "a") == "x"

        # line endings are normalized as XML 1.0 §2.11 requires
        @test MortalityTables._content(node("<a>a\r\nb</a>")) == "a\nb"
        @test MortalityTables._content(node("<a>x &amp; y</a>")) == "x & y"
    end

    @testset "XTbML loading" begin
        pth = joinpath(soa_tbl_dir,"t1076.xml")
        xtbl = MortalityTables.parseXTbMLTable(MortalityTables.open_and_read(pth), pth)
        @test xtbl isa NamedTuple
        @test xtbl.metadata.id == "1076"
        @test xtbl.metadata.content_type == "CSO/CET"
    end

    @testset "readXTbML cache" begin
        pth = joinpath(soa_tbl_dir,"t1076.xml")
        # repeated reads of the same path return the identical object
        @test MortalityTables.readXTbML(pth) === MortalityTables.readXTbML(pth)
    end

    @testset "Ultimate Only" begin
        mt = MortalityTables._read_xtbml(joinpath(soa_tbl_dir,"t17.xml"))
        @test isa(mt, MortalityTable)

        @test mt.ultimate[0] ≈ 0.00245
        @test mt.ultimate[100] ≈ 1.0
        @test_throws BoundsError mt.ultimate[101]
    end

    @testset "XTbML to MortalityTable" begin
        @testset "Select and Ultimate" begin
            mt = MortalityTables._read_xtbml(joinpath(soa_tbl_dir,"t1076.xml"))
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
            build(sel, u = ult) = MortalityTables._table_from_labels((select = sel, ultimate = u, metadata = md), "XTbML table probe")
            # consecutive durations run into the ultimate rates at the next attained age
            mt = build(row(1:3))
            @test mt.select[40][40:44] == [0.01, 0.02, 0.03, 0.2, 0.2]
            @test eltype(mt.select[40]) == Float64
            # a missing leading or interior duration keeps every other rate at its own age
            @test isequal(build(row(2:3)).select[40][40:43], [missing, 0.02, 0.03, 0.2])
            @test isequal(build(row([1, 3])).select[40][40:43], [0.01, missing, 0.03, 0.2])
            # an empty cell is a missing rate: an interior one stays missing, and trailing ones end
            # the select period, as when the durations are absent
            cells = [(duration = d, rate = r) for (d, r) in zip(1:5, [0.01, missing, 0.03, missing, missing])]
            @test isequal(build([(issue_age = 40, rates = cells)]).select[40], build(row([1, 3])).select[40])
            @test eltype(build([(issue_age = 40, rates = cells[1:1])]).select[40]) == Float64
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
            mt = MortalityTables._read_xtbml(joinpath(soa_tbl_dir,"t18.xml"))
            @test isa(mt, MortalityTable)
        end

        @testset "empty metadata elements" begin
            # t217 has an empty metadata element that the 2.x parser could not read
            @test MortalityTables.table(217) isa MortalityTables.UltimateTable
        end

        @testset "layouts with more than two tables are rejected" begin
            # t357 has three <Table> elements: reading its first and last would silently drop
            # the issue ages in the middle one
            path = joinpath(soa_tbl_dir, "t357.xml")
            @test_throws ArgumentError MortalityTables._read_xtbml(path)
            @test_throws "has 3 <Table> elements" MortalityTables._read_xtbml(path)
        end

        @testset "every bundled table loads, except the known unsupported ones" begin
            # test/data/unsupported_tables.txt lists each table that does not load, with its error
            expected = Dict{Int, String}()
            for line in eachline(joinpath(pkgdir(MortalityTables), "test", "data", "unsupported_tables.txt"))
                (isempty(line) || startswith(line, "#")) && continue
                id, err = split(line)
                expected[parse(Int, id)] = err
            end
            files = filter(f -> endswith(f, ".xml") && !startswith(f, "._"), readdir(soa_tbl_dir))
            failed = Dict{Int, String}()
            for f in files
                try
                    # uncached, so this does not fill the table cache
                    MortalityTables._read_xtbml(joinpath(soa_tbl_dir, f))
                catch e
                    failed[parse(Int, f[2:(end - 4)])] = string(nameof(typeof(e)))
                end
            end
            @test length(files) - length(failed) > 2000
            # the ids of tables that newly fail or newly load
            @test sort!(collect(setdiff(keys(failed), keys(expected)))) == Int[]
            @test sort!(collect(setdiff(keys(expected), keys(failed)))) == Int[]
            @test [id for id in keys(failed) if get(expected, id, nothing) != failed[id]] == Int[]
        end
    end
end
