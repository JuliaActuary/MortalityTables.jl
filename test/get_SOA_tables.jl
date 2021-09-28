@testset "HTTP requests to mort.SOA.org" begin

    tbl = get_SOA_table(60029)

    @test tbl isa MortalityTable
    @test tbl[0] == 0.10139


    get_SOA_table(60029)

    @test MortalityTables.table("Australian Life Tables 1891-1900 Female") isa MortalityTable
    @test MortalityTables.table(60029) isa MortalityTable
    @test MortalityTables.table("Australian Life Tables 1891-1900 Female")[0] == 0.10139
    
    # issue #69
    tbl = get_SOA_table(887)
    @test tbl isa MortalityTable 
    tbl = get_SOA_table(2585)
    @test tbl isa MortalityTable 

    @test_throws ArgumentError get_SOA_table("hello")
end
