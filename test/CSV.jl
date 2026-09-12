tbl_dir_test = joinpath(pkgdir(MortalityTables), "test", "data", "CSV")

@testset "CSV select row parsing" begin
    ext = Base.get_extension(MortalityTables, :MortalityTablesCSVExt)
    @test ext !== nothing
    # trailing blanks are trimmed; an interior blank is kept, not truncated at
    rates = ext._select_rates([0.1, missing, 0.3, missing, missing])
    @test length(rates) == 3
    @test rates[1] == 0.1
    @test rates[2] === missing
    @test rates[3] == 0.3
    # string cells are parsed, and a row with no blanks stays Float64
    @test ext._select_rates(["0.1", "0.2"]) == [0.1, 0.2]
    @test eltype(ext._select_rates([0.1, 0.2, missing])) == Float64
    @test isempty(ext._select_rates([missing, missing]))
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