tbl_dir = joinpath(pkgdir(MortalityTables),"src","tables","SOA")

@testset "CSV and XTbML equality: $id" for id in [17,428,1152,3302]
    xtbml = MortalityTables.readXTbML(joinpath(tbl_dir,"t$id.xml"))
    csv = MortalityTables.readcsv(joinpath(tbl_dir,"t$id.csv"))
    
    @test xtbml.ultimate == csv.ultimate
    if typeof(xtbml) <: MortalityTables.SelectUltimateTable
        for issue_age in eachindex(xtbml.select)
            @test xtbml.select[issue_age] == csv.select[issue_age]
        end
    end

end