tbl_dir_test = joinpath(pkgdir(MortalityTables), "test", "data", "CSV")
@testset "CSV equality" begin
    @testset "CSV and XTbML equality: $id" for id in [17,428,1152,3302]
        xtbml = MortalityTables.readXTbML(joinpath(tbl_dir,"SOA", "t$id.xml"))
        csv = MortalityTable(CSV.File(joinpath(tbl_dir_test, "t$id.csv"), header=false, silencewarnings=true))
        
        @test xtbml.ultimate == csv.ultimate
        if typeof(xtbml) <: MortalityTables.SelectUltimateTable
            for issue_age in eachindex(xtbml.select)
                @test xtbml.select[issue_age] == csv.select[issue_age]
            end
        end

    end
end