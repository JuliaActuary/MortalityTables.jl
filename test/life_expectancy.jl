@testset "life expectancy" begin
    
    @testset "parametric" begin

        # ALMC 2.5
        m = MortalityTables.Gompertz(a=0.0003,b=log(1.07))

        @test life_expectancy(m,0) ≈ 71.938 atol=0.001
        @test life_expectancy(m,50) ≈ 26.691 atol=0.001

        # a DeathDistribution is accepted and ignored by continuous models
        @test life_expectancy(m, 50, MortalityTables.UniformDeaths()) == life_expectancy(m, 50)

        # a law defined at every age (omega = Inf) gives the same bits as integrating the
        # remaining lifetime from zero, as before omega bounded the integral
        for law in (m, Makeham()), age in (0, 40, 50.5)
            @test life_expectancy(law, age) === quadgk(to -> survival(law, age, to + age), 0, Inf)[1]
        end
    end

    @testset "parametric, to omega" begin
        # Wittstein's formula ends at m = 100, where survival from 40 is still about 4.7e-6.
        # Reference: 256-bit quadrature of the hazard written out independently.
        w = MortalityTables.Wittstein()
        @test life_expectancy(w, 40) ≈ 8.480415927762158 rtol = 1.0e-8

        # a hazard of 1/(100 - age) makes the remaining lifetime from 40 uniform on [40, 100]
        v = MortalityTables.VanderMaen(a = 0.0, b = 0.0, c = 0.0, i = 1.0, n = 100.0)
        v2 = MortalityTables.VanderMaen2(a = 0.0, b = 0.0, i = 1.0, n = 100.0)
        @test life_expectancy(v, 40) ≈ 30
        @test life_expectancy(v2, 40) ≈ 30

        # nothing remains at omega, in the type of the integral
        for law in (w, v, v2), age in (100, 100.0)
            @test life_expectancy(law, age) === 0.0
        end
        w32 = MortalityTables.Wittstein(a = 1.5f0, b = 1.0f0, n = 0.5f0, m = 100.0f0)
        @test life_expectancy(w32, 40.0f0) isa Float32
        @test life_expectancy(w32, 100.0f0) === 0.0f0

        # past omega the law ends, and life_expectancy throws
        @test_throws DomainError life_expectancy(w, 101)
        @test_throws DomainError life_expectancy(v, 101)
        @test_throws DomainError life_expectancy(v2, 101)
        # also where Wittstein's formula still evaluates (a whole-number n), and just past omega
        for n in (1.0, 2.0)
            wn = MortalityTables.Wittstein(n = n)
            for age in (101, 100.01, 100 + 1.0e-9)
                @test_throws DomainError life_expectancy(wn, age)
            end
            # at and below omega, the integral of survival as before
            H(x) = cumhazard(wn, x)
            @test life_expectancy(wn, 99) === quadgk(t -> exp(-(H(99 + t) - H(99))), 0, omega(wn) - 99)[1]
            @test life_expectancy(wn, 100) === 0.0
        end
    end

    @testset "vector" begin
        # if no distribution of deaths given, assume curtate otherwise complete
        t = MortalityTables.table("1980 CSO Basic Table – Male, ANB")

        # calculate sum of tpx in Excel
        @test life_expectancy(t.ultimate,55) ≈ 22.16469212 atol=1e-6
        @test life_expectancy(t.ultimate,100) ≈ 0.0 atol=1e-6
        # the table has no rates past omega
        @test_throws BoundsError life_expectancy(t.ultimate,101)
        @test_throws BoundsError life_expectancy(t.ultimate,-1)
        # past omega, the integral runs back through ages the table does not have
        @test_throws BoundsError life_expectancy(t.ultimate,102,MortalityTables.UniformDeaths())
        # values at several ages match the quadratic formulation this replaced
        for age in (0, 35, 55, 99)
            @test life_expectancy(t.ultimate, age) ≈ sum(survival(t.ultimate, age, age + dur) for dur in 1:omega(t.ultimate)-age)
        end

        @test_throws MethodError life_expectancy(t,100)
        @test_throws MethodError life_expectancy(t,100,MortalityTables.UniformDeaths())

        # relation of curtate to complete, ALMC 2.6.1
        @test life_expectancy(t.ultimate,55,MortalityTables.UniformDeaths()) ≈ 22.16469212 + 0.5 atol=1e-3
    end

end