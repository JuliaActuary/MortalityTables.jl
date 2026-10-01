@testset "life expectancy" begin
    
    @testset "parametric" begin

        # ALMC 2.5
        m = MortalityTables.Gompertz(a=0.0003,b=log(1.07))

        @test life_expectancy(m,0) ≈ 71.938 atol=0.001
        @test life_expectancy(m,50) ≈ 26.691 atol=0.001

        # a DeathDistribution is accepted and ignored by continuous models
        @test life_expectancy(m, 50, MortalityTables.Uniform()) == life_expectancy(m, 50)
    end

    @testset "vector" begin
        # if no distribution of deaths given, assume curtate otherwise complete
        t = MortalityTables.table("1980 CSO Basic Table – Male, ANB")

        # calculate sum of tpx in Excel
        @test life_expectancy(t.ultimate,55) ≈ 22.16469212 atol=1e-6
        @test life_expectancy(t.ultimate,100) ≈ 0.0 atol=1e-6
        # an age past omega has no remaining years
        @test life_expectancy(t.ultimate,101) == 0.0
        # values at several ages match the quadratic formulation this replaced
        for age in (0, 35, 55, 99)
            @test life_expectancy(t.ultimate, age) ≈ sum(survival(t.ultimate, age, age + dur) for dur in 1:omega(t.ultimate)-age)
        end

        @test_throws ArgumentError life_expectancy(t,100)

        # relation of curtate to complete, ALMC 2.6.1
        @test life_expectancy(t.ultimate,55,MortalityTables.Uniform()) ≈ 22.16469212 + 0.5 atol=1e-3
    end

end