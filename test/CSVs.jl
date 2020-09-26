@testset "CSV and XTbML equality" begin
    tbl_dir = joinpath(pkgdir(MortalityTables),"src","tables","SOA")

    for id in [17,428,1152,3302]
        xtbml = MortalityTables.readXTbML(joinpath(tbl_dir,"t$id.xml"))
        csv = MortalityTables.readcsv(joinpath(tbl_dir,"t$id.csv"))

        @show id
        @show xtbml.ultimate, firstindex(xtbml.ultimate)
        @show csv.ultimate, firstindex(csv.ultimate)
        @test xtbml.ultimate == csv.ultimate
    end


end