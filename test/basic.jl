
@testset "basic `Ultimate` MortalityTable" begin

    q1 = UltimateMortality([i for i = 1:20], 0) |> rates
    @test q1[0] == 1
    @test q1[1] == 2
    @test q1[0:1] == [1, 2]
    @test ω(q1) == 19
    @test_throws BoundsError q1[ω(q1) + 1]

    # non-zero start age
    q2 = UltimateMortality([i for i = 1:10], 5) |> rates
    @test_throws BoundsError q2[4]
    @test q2[5] == 1

    # select strucutre
    s = [ia + d for ia = 0:5, d = 1:5]

    s1 = SelectMortality(s, u1, 0)
    q3 = rates(s1,0)
    @test q3[0] == 1
    @test q3[0:1] == [1, 2]
    @test ω(q3) == 19
    @test_throws BoundsError q3[ω(q3) + 1]
    @test q3[19] == 20
    
    q4 = rates(s1,5)
    @test q4[6] == 5


    mt1 = MortalityTable(s1, u1, TableMetaData())

    @test mt1.select[0] == 1
    @test mt1.ultimate[0][1] == 1

    # test time zero accumlated force
    @test cumulative_decrement(mt1.ultimate[0],0) == 0
    @test survivorship(mt1.ultimate[0],0) == 1

end
