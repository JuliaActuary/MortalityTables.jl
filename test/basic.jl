
@testset "basic MortalityTable" begin

    q1 = UltimateMortality([i for i = 0:19], 0)
    @test q1[0] == 0
    @test q1[1] == 1
    @test q1[0:1] == [0, 1]
    @test omega(q1) == 19
    @test_throws BoundsError q1[omega(q1) + 1]

    # non-zero start age
    q2 = UltimateMortality([i for i = 0:9], 5)
    @test_throws BoundsError q2[4]
    @test q2[5] == 0

    # select strucutre
    s = [ia + d for ia = 0:5, d = 0:4]

    q3 = SelectMortality(s, q1, 0)
    @test q3[0][0] == 0
    @test q3[0][0:1] == [0, 1]
    @test omega(q3[0]) == 19
    @test_throws BoundsError q3[omega(q3[0]) + 1]
    @test q3[0][19] == 19
    
    @test q3[5][5] == 5


    mt1 = MortalityTable(q3, q1, TableMetaData())

    @test mt1.select[0][1] == 1
    @test mt1.ultimate[1] == 1

    mt2 = MortalityTable(q1,TableMetaData())

    @test mt2.ultimate[0] == 0
    @test mt2[0] == 0

    # test time zero accumlated force

    q4 = UltimateMortality([0.1,0.3,0.6,1],0)
    
    @test survivorship(q4,0) ≈ 1
    @test cumulative_decrement(q4,0) ≈ 0

    @test survivorship(q4,1) ≈ 0.9
    @test cumulative_decrement(q4,1) ≈ 0.1

    @test survivorship(q4,1,1) ≈ 1.0
    @test survivorship(q4,1,2) ≈ 0.7
    @test cumulative_decrement(q4,1,1) ≈ 0.0
    @test cumulative_decrement(q4,1,2) ≈ 0.3
    
    @test survivorship(q4,1,4) ≈ 0.0
    @test cumulative_decrement(q4,1,4) ≈ 1.0

end
