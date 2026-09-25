tbl_dir_test = joinpath(pkgdir(MortalityTables), "test", "data", "CSV")

@testset "CSV select row parsing" begin
    ext = Base.get_extension(MortalityTables, :MortalityTablesCSVExt)
    @test ext !== nothing
    # trailing blanks are trimmed; an interior blank is kept, not truncated at
    rates = ext._select_rates([0.1, missing, 0.3, missing, missing], 1:5)
    @test length(rates) == 3
    @test rates[1] == 0.1
    @test rates[2] === missing
    @test rates[3] == 0.3
    # string cells are parsed, and a row with no blanks stays Float64
    @test ext._select_rates(["0.1", "0.2"], 1:2) == [0.1, 0.2]
    @test eltype(ext._select_rates([0.1, 0.2, missing], 1:3)) == Float64
    @test isempty(ext._select_rates([missing, missing], 1:2))
    # rates go to their duration labels, not their column positions
    @test isequal(ext._select_rates([0.1, 0.2, 0.4], [1, 2, 4]), [0.1, 0.2, missing, 0.4])
    @test ext._select_rates([0.2, 0.1], [2, 1]) == [0.1, 0.2]
    @test_throws "repeat" ext._select_rates([0.1, 0.2], [1, 1])
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

    # the same labels read as XTbML give the same table
    md = MortalityTables.TableMetaData(name = "gapped select")
    xml = MortalityTables.XTbML_Table_To_MortalityTable(
        (
            select = [
                (issue_age = 40, rates = [(duration = d, rate = r) for (d, r) in zip([1, 2, 4], [0.01, 0.02, 0.04])]),
                (issue_age = 42, rates = [(duration = d, rate = r) for (d, r) in zip([1, 4], [0.03, 0.05])]),
            ],
            ultimate = [(age = a, rate = r) for (a, r) in zip(40:46, [0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26])],
            metadata = md,
        )
    )
    @test isequal(csv.ultimate, xml.ultimate)
    @test axes(csv.select) == axes(xml.select)
    for issue_age in eachindex(xml.select)
        @test isequal(csv.select[issue_age], xml.select[issue_age])
    end

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