
@testset "basic MortalityTable" begin
    @testset "basic structure" begin
        q1 = UltimateMortality([i for i = 0:19], start_age = 0)
        @test q1[0] == 0
        @test q1[1] == 1
        @test q1[0:1] == [0, 1]
        @test omega(q1) == 19
        @test_throws BoundsError q1[omega(q1) + 1]

        # non-zero start age
        q2 = UltimateMortality([i for i = 0:9], start_age = 5)
        @test_throws BoundsError q2[4]
        @test q2[5] == 0

        # select strucutre
        s = [ia + d for ia = 0:5, d = 0:4]

        q3 = SelectMortality(s, q1, start_age = 0)
        @test q3[0][0] == 0
        @test q3[0][0:1] == [0, 1]
        @test omega(q3[0]) == 19
        @test_throws BoundsError q3[omega(q3[0]) + 1]
        @test q3[0][19] == 19
        
        @test q3[5][5] == 5


        mt1 = MortalityTable(q3, q1)

        @test mt1.select[0][1] == 1
        @test mt1.ultimate[1] == 1

        mt2 = MortalityTable(q1)

        @test mt2.ultimate[0] == 0
        @test mt2[0] == 0

    end

    @testset "off-aligned select and ult" begin
        select_matrix = [(i + j - 1) / 100 for i in 0:10,j in 1:20]
        ult = UltimateMortality([i / 100 for i in 18:100],start_age=18)

        select = SelectMortality(select_matrix,ult)

        for issue_age in 0:10
            for dur in 1:30
                @test select[issue_age][issue_age + dur - 1] == (issue_age + dur - 1) / 100
            end
        end
    end

    # test time zero accumlated force
    @testset "accumulated force" begin
        q4 = UltimateMortality([0.1,0.3,0.6,1])
        
        @test survival(q4, 0) ≈ 1
        @test decrement(q4, 0) ≈ 0

        @test survival(q4, 1) ≈ 0.9
        @test decrement(q4, 1) ≈ 0.1

        @test survival(q4, 1, 1) ≈ 1.0
        @test survival(q4, 1, 2) ≈ 0.7
        @test decrement(q4, 1, 1) ≈ 0.0
        @test decrement(q4, 1, 2) ≈ 0.3
        
        @test survival(q4, 1, 4) ≈ 0.0
        @test decrement(q4, 1, 4) ≈ 1.0

        
        @test survival(q4, -1) ≈ 1.0
        @test survival(q4, 4, -1) ≈ 1.0
        @test decrement(q4, -1) ≈ 0.0
        @test decrement(q4, 4, -1) ≈ 0.0
    end

    @testset "Metadata" begin
        d = TableMetaData()

        @test isnothing(d.name)
        
        d = TableMetaData(name = "test")
        @test d.name == "test"
    end

    @testset "mortality_vector" begin
        v = [i for i in 3:10]
        q = mortality_vector(v,start_age=3)
        @test q[3] == 3
        @test q[10] == 10

        q = mortality_vector(collect(0:5))
        @test q[0] == 0
        @test q[5] == 5
    end

end
