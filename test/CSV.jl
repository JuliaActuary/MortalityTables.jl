tbl_dir_test = joinpath(pkgdir(MortalityTables), "test", "data", "CSV")

@testset "CSV select row records" begin
    ext = Base.get_extension(MortalityTables, :MortalityTablesCSVExt)
    @test ext !== nothing
    # a blank cell is a missing rate, string cells are parsed, and each rate keeps its duration label
    rec(cells, durations) = ext._select_records(cells, durations)
    @test isequal(
        rec([0.1, missing, 0.3, missing, missing], 1:5),
        [(duration = d, rate = r) for (d, r) in zip(1:5, [0.1, missing, 0.3, missing, missing])]
    )
    @test rec(["0.1", "0.2"], 1:2) == [(duration = 1, rate = 0.1), (duration = 2, rate = 0.2)]
    @test rec([0.1, 0.2, 0.4], [1, 2, 4]) == [(duration = 1, rate = 0.1), (duration = 2, rate = 0.2), (duration = 4, rate = 0.4)]
    # placed by the shared assembler: trailing blanks shorten the row and an interior blank is kept
    md = MortalityTables.TableMetaData(name = "records")
    ult = [(age = a, rate = 0.5) for a in 40:50]
    build(rates) = MortalityTables._table_from_labels((select = [(issue_age = 40, rates)], ultimate = ult, metadata = md), "CSV table records")
    @test isequal(build(rec([0.1, missing, 0.3, missing, missing], 1:5)).select[40][40:43], [0.1, missing, 0.3, 0.5])
    @test eltype(build(rec([0.1, 0.2, missing], 1:3)).select[40]) == Float64
    @test build(rec([0.2, 0.1], [2, 1])).select[40][40:41] == [0.1, 0.2]
    @test_throws "repeat" build(rec([0.1, 0.2], [1, 1]))
end

@testset "CSV rates are placed by their labels" begin
    read_csv(text) = MortalityTable(CSV.File(IOBuffer(text); header = false, silencewarnings = true))

    # ultimate ages 60 and 62: the rate labelled 62 stays at 62
    mt = read_csv("""
    Table Name:,gapped ultimate
    Table Identity:,999

    Row\\Column,Rate
    60,0.1
    62,0.3
    """)
    @test axes(mt.ultimate, 1) == 60:62
    @test mt.ultimate[60] == 0.1 && mt.ultimate[62] == 0.3
    @test ismissing(mt.ultimate[61])

    # grouped issue ages (40, 42), nonconsecutive duration headers (1, 2, 4), and a blank cell
    csv = read_csv("""
    Table Name:,gapped select,,
    Table Identity:,998,,

    Row\\Column,1,2,4
    40,0.01,0.02,0.04
    42,0.03,,0.05

    Table # ,2,,
    Row\\Column,1,,
    40,0.2,,
    41,0.21,,
    42,0.22,,
    43,0.23,,
    44,0.24,,
    45,0.25,,
    46,0.26,,
    """)
    @test isequal(csv.select[40][40:44], [0.01, 0.02, missing, 0.04, 0.24])
    @test isequal(csv.select[42][42:46], [0.03, missing, missing, 0.05, 0.26])
    @test ismissing(csv.select[41])

    # a select row with no rates is an absent issue age, as in XTbML, not a zero-length select
    # period that falls through to the ultimate rates
    blank = read_csv("""
    Table Name:,blank select row,,
    Table Identity:,996,,

    Row\\Column,1,2,4
    40,0.01,0.02,0.04
    41,,,
    42,0.03,,0.05

    Table # ,2,,
    Row\\Column,1,,
    40,0.2,,
    41,0.21,,
    42,0.22,,
    43,0.23,,
    44,0.24,,
    45,0.25,,
    46,0.26,,
    """)
    @test ismissing(blank.select[41])
    @test axes(blank.select) == axes(csv.select)
    for issue_age in (40, 42)
        @test isequal(blank.select[issue_age], csv.select[issue_age])
    end

    # the same labels, as XTbML parses them, give the same table
    md = MortalityTables.TableMetaData(name = "gapped select")
    xml = MortalityTables._table_from_labels(
        (
            select = [
                (issue_age = 40, rates = [(duration = d, rate = r) for (d, r) in zip([1, 2, 4], [0.01, 0.02, 0.04])]),
                (issue_age = 42, rates = [(duration = d, rate = r) for (d, r) in zip([1, 4], [0.03, 0.05])]),
            ],
            ultimate = [(age = a, rate = r) for (a, r) in zip(40:46, [0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26])],
            metadata = md,
        ),
        "XTbML table gapped select"
    )
    @test isequal(csv.ultimate, xml.ultimate)
    @test axes(csv.select) == axes(xml.select)
    for issue_age in eachindex(xml.select)
        @test isequal(csv.select[issue_age], xml.select[issue_age])
    end

    # a select period that runs past the ultimate omega has no ultimate tail, read either way
    past = read_csv("""
    Table Name:,long select,,
    Table Identity:,995,,

    Row\\Column,1,2,3
    44,0.01,0.02,0.03
    45,0.04,0.05,0.06

    Table # ,2,,
    Row\\Column,1,,
    44,0.2,,
    45,0.21,,
    46,0.22,,
    """)
    past_xml = MortalityTables._table_from_labels(
        (
            select = [
                (issue_age = a, rates = [(duration = d, rate = r) for (d, r) in zip(1:3, rs)])
                    for (a, rs) in ((44, [0.01, 0.02, 0.03]), (45, [0.04, 0.05, 0.06]))
            ],
            ultimate = [(age = a, rate = r) for (a, r) in zip(44:46, [0.2, 0.21, 0.22])],
            metadata = MortalityTables.TableMetaData(name = "long select"),
        ),
        "XTbML table long select"
    )
    @test axes(past.select[45], 1) == 45:47
    @test past.select[45][47] == 0.06
    @test isequal(past.ultimate, past_xml.ultimate)
    for issue_age in (44, 45)
        @test isequal(past.select[issue_age], past_xml.select[issue_age])
    end

    # a third block of rates would be dropped
    @test_throws "has 3 tables" read_csv("""
    Table Name:,three tables

    Row\\Column,Rate
    60,0.1

    Row\\Column,Rate
    60,0.2

    Row\\Column,Rate
    60,0.3
    """)

    # a repeated age label is ambiguous
    @test_throws "repeat" read_csv("""
    Table Name:,repeated ultimate
    Table Identity:,997

    Row\\Column,Rate
    60,0.1
    60,0.3
    """)
end

@testset "CSV equality" begin
    @testset "CSV and XTbML equality: $id" for id in [17,428,1152,3302]
        xtbml = MortalityTables.readXTbML(joinpath(soa_tbl_dir, "t$id.xml"))
        csv = MortalityTable(CSV.File(joinpath(tbl_dir_test, "t$id.csv"), header=false, silencewarnings=true))
        
        @test xtbml.ultimate == csv.ultimate
        if typeof(xtbml) <: MortalityTables.SelectUltimateTable
            for issue_age in eachindex(xtbml.select)
                @test xtbml.select[issue_age] == csv.select[issue_age]
            end
        end

    end
end