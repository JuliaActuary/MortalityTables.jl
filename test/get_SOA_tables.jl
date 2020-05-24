@testset "HTTP requests to mort.SOA.org" begin
    tbl = get_SOA_table(60029)

    @test tbl isa MortalityTable
    @test q(tbl,0,1) == 0.10139


    get_SOA_table!(tables,60029)

    @test tables["Australian Life Tables 1891-1900 Female"] isa MortalityTable
    @test q(tables["Australian Life Tables 1891-1900 Female"],0,1) == 0.10139
end
