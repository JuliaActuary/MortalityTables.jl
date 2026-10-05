# A law with survival S(t) = (1 + t)^-p from age 0: hazard p / (1 + age). For p = 2 the curtate
# expectancy at age 0 is Σ_{k≥1} (1 + k)^-2 = π²/6 - 1, and its terms fall off slowly.
struct InversePower{T} <: MortalityTables.ParametricMortality
    p::T
end
MortalityTables.hazard(m::InversePower, age) = m.p / (1 + age)
MortalityTables.cumhazard(m::InversePower, age) = m.p * log1p(age)

@testset "life expectancy" begin

    @testset "parametric, omega = Inf" begin
        # ALMC 2.5
        m = MortalityTables.Gompertz(a = 0.0003, b = log(1.07))

        @test complete_life_expectancy(m, 0) ≈ 71.938 atol = 0.001
        @test complete_life_expectancy(m, 50) ≈ 26.691 atol = 0.001

        # the curtate sum, against its definition summed in BigFloat until the terms vanish
        mbig = MortalityTables.Gompertz(a = big"0.0003", b = log(big"1.07"))
        for age in (0, 50, 50.5)
            ref = sum(survival(mbig, big(age), big(age) + k) for k in 1:400)
            @test curtate_life_expectancy(m, age) ≈ ref rtol = sqrt(eps())
            @test (@inferred curtate_life_expectancy(m, age)) isa Float64
        end

        # a DeathDistribution is accepted and ignored by continuous models
        @test complete_life_expectancy(m, 50, UniformDeaths()) === complete_life_expectancy(m, 50)
        @test complete_life_expectancy(m, 50, ConstantForce()) === complete_life_expectancy(m, 50)

        # a law defined at every age (omega = Inf) gives the same bits as integrating the
        # remaining lifetime from zero, as before omega bounded the integral
        for law in (m, Makeham()), age in (0, 40, 50.5)
            @test complete_life_expectancy(law, age) === quadgk(to -> survival(law, age, to + age), 0, Inf)[1]
        end
    end

    @testset "parametric curtate tail, omega = Inf" begin
        # the terms of Σ (1 + k)^-2 fall off slowly: a rule that stops once the next term is small
        # gave 0.6446562. The midpoint of the tail's bracket meets the tolerance.
        exact = π^2 / 6 - 1
        for T in (Float64, Float32)
            e = @inferred curtate_life_expectancy(InversePower(T(2)), zero(T))
            @test e isa T
            @test abs(e - exact) <= sqrt(eps(T)) * exact
        end
        # the tolerances are keywords
        m = InversePower(2.0)
        @test abs(curtate_life_expectancy(m, 0.0; rtol = 0, atol = 1.0e-3) - exact) <= 1.0e-3
        @test abs(curtate_life_expectancy(m, 0.0; rtol = 1.0e-10) - exact) <= 1.0e-10 * exact
        # a tolerance out of reach within 100,000 terms does not converge, rather than truncate
        @test_throws ArgumentError curtate_life_expectancy(m, 0.0; rtol = 1.0e-12)
        # so does a law whose survival levels off (the hazard tends to zero), whose curtate
        # expectancy is infinite
        @test_throws ArgumentError curtate_life_expectancy(InversePower(0.5), 0.0)
        # from a later age the sum is conditional on survival to that age
        @test curtate_life_expectancy(m, 1.0) ≈ 4 * (π^2 / 6 - 1 - 1 / 4) rtol = sqrt(eps())
    end

    @testset "parametric complete expectancy converges or throws" begin
        # S(t) = (1 + t)^-p integrates to 1 / (p - 1) for p > 1; from age 1, (2 / (2 + t))^3 gives 1
        @test complete_life_expectancy(InversePower(3.0), 0.0) ≈ 0.5 rtol = sqrt(eps())
        @test complete_life_expectancy(InversePower(3.0), 1.0) ≈ 1.0 rtol = sqrt(eps())
        @test abs(complete_life_expectancy(InversePower(3.0), 0.0; rtol = 1.0e-10) - 0.5) <= 1.0e-10
        # for p ≤ 1 the integral is infinite, and QuadGK's estimate does not converge
        @test_throws ArgumentError complete_life_expectancy(InversePower(1.0), 0.0)
        @test_throws ArgumentError complete_life_expectancy(InversePower(0.5), 0.0)
        # InverseWeibull's survival also falls too slowly: QuadGK returned about 6.5e8 years with an
        # error estimate of 1.7e8
        @test_throws ArgumentError complete_life_expectancy(MortalityTables.InverseWeibull(), 30)
        # a finite integral with too few evaluations does not converge either
        m = MortalityTables.Gompertz(a = 0.0003, b = log(1.07))
        @test_throws ArgumentError complete_life_expectancy(m, 50; maxevals = 1)
        @test complete_life_expectancy(m, 50; maxevals = 10^7) === complete_life_expectancy(m, 50)
        # the curtate sum needs finite estimates too
        @test_throws ArgumentError curtate_life_expectancy(InversePower(1.0), 0.0)
    end

    @testset "parametric, to a finite omega" begin
        # Wittstein's formula ends at m = 100, where survival from 40 is still about 4.7e-6.
        # Reference: 256-bit quadrature of the hazard written out independently.
        w = MortalityTables.Wittstein()
        @test complete_life_expectancy(w, 40) ≈ 8.480415927762158 rtol = 1.0e-8

        # a hazard of 1/(100 - age) makes the remaining lifetime from 40 uniform on [40, 100]
        v = MortalityTables.VanderMaen(a = 0.0, b = 0.0, c = 0.0, i = 1.0, n = 100.0)
        v2 = MortalityTables.VanderMaen2(a = 0.0, b = 0.0, i = 1.0, n = 100.0)
        for law in (v, v2)
            @test complete_life_expectancy(law, 40) ≈ 30
            # omega - age is whole: the last term is the survival to omega, which here is zero.
            # Σ_{k=1}^{60} (60 - k)/60
            @test curtate_life_expectancy(law, 40) ≈ 29.5
            # omega - age is fractional: the sum stops at ⌊omega - age⌋ = 59.
            # Σ_{k=1}^{59} (59.5 - k)/59.5
            @test curtate_life_expectancy(law, 40.5) ≈ 59 - 1770 / 59.5
        end
        # survival left at omega counts as death there: the last whole year to omega counts in
        # the curtate sum, and a fraction of a year before omega gives no whole year
        @test survival(w, 99, 100) > 0
        @test curtate_life_expectancy(w, 99) === survival(w, 99, 100)
        @test curtate_life_expectancy(w, 98.5) === survival(w, 98.5, 99.5)
        @test curtate_life_expectancy(w, 99.5) === 0.0
        @test curtate_life_expectancy(w, 40) ≈ sum(survival(w, 40, 40 + k) for k in 1:60)

        # nothing remains at omega, in the type of the result
        for law in (w, v, v2), age in (100, 100.0)
            @test complete_life_expectancy(law, age) === 0.0
            @test curtate_life_expectancy(law, age) === 0.0
        end
        w32 = MortalityTables.Wittstein(a = 1.5f0, b = 1.0f0, n = 0.5f0, m = 100.0f0)
        @test complete_life_expectancy(w32, 40.0f0) isa Float32
        @test complete_life_expectancy(w32, 100.0f0) === 0.0f0
        @test curtate_life_expectancy(w32, 40.0f0) isa Float32
        @test curtate_life_expectancy(w32, 100.0f0) === 0.0f0

        # past omega the law ends, and both throw, also just past omega where the curtate sum
        # would be empty and the integral short
        for law in (w, v, v2), age in (101, 100.5, 100 + 1.0e-9)
            @test_throws DomainError complete_life_expectancy(law, age)
            @test_throws DomainError curtate_life_expectancy(law, age)
            @test_throws DomainError curtate_life_expectancy(law, age; rtol = 1.0e-3)
        end
        # also where Wittstein's formula still evaluates (a whole-number n)
        for n in (1.0, 2.0)
            wn = MortalityTables.Wittstein(n = n)
            for age in (101, 100.01, 100 + 1.0e-9)
                @test_throws DomainError complete_life_expectancy(wn, age)
                @test_throws DomainError curtate_life_expectancy(wn, age)
            end
            # at and below omega, the integral of survival as before
            H(x) = cumhazard(wn, x)
            @test complete_life_expectancy(wn, 99) === quadgk(t -> exp(-(H(99 + t) - H(99))), 0, omega(wn) - 99)[1]
            @test complete_life_expectancy(wn, 100) === 0.0
        end
    end

    @testset "vector" begin
        t = MortalityTables.table("1980 CSO Basic Table – Male, ANB")
        q = t.ultimate

        # calculate sum of tpx in Excel
        @test curtate_life_expectancy(q, 55) ≈ 22.16469212 atol = 1.0e-6
        @test curtate_life_expectancy(q, 100) === 0.0
        @test complete_life_expectancy(q, 100) === 0.0
        # an age outside the table is a BoundsError, also where the sum or the integral
        # would otherwise be empty or run through a year the table has
        for age in (101, 102, -1)
            @test_throws BoundsError curtate_life_expectancy(q, age)
            @test_throws BoundsError complete_life_expectancy(q, age)
        end
        for age in (100.5, 100 + 1.0e-9, -0.5), dd in (UniformDeaths(), ConstantForce(), Balducci())
            @test_throws BoundsError complete_life_expectancy(q, age, dd)
        end
        # values at several ages match the quadratic formulation this replaced
        for age in (0, 35, 55, 99)
            @test curtate_life_expectancy(q, age) ≈ sum(survival(q, age, age + dur) for dur in 1:(omega(q) - age))
        end

        # a whole table is not a vector of rates
        @test_throws MethodError curtate_life_expectancy(t, 100)
        @test_throws MethodError complete_life_expectancy(t, 100)
        @test_throws MethodError complete_life_expectancy(t, 100, UniformDeaths())

        # relation of curtate to complete, ALMC 2.6.1
        @test complete_life_expectancy(q, 55) ≈ 22.16469212 + 0.5 atol = 1.0e-3
        # under uniform deaths it is exact: the expectancy to the last age is the curtate one
        # plus half of the probability of dying before the last age
        for age in (0, 35, 55, 99, 100)
            @test complete_life_expectancy(q, age, UniformDeaths()) ≈
                curtate_life_expectancy(q, age) + (1 - survival(q, age, omega(q))) / 2 rtol = 1.0e-14
        end
    end

    @testset "vector, fractional ages" begin
        qs = UltimateMortality([0.1, 0.3, 1.0])
        # by hand under uniform deaths: from 0.5 to 1 survival is (1 - u·0.1)/0.95, with integral
        # 0.4625/0.95, then 0.9/0.95 survive to 1 and live 1 - 0.3/2 of the next year
        @test complete_life_expectancy(qs, 0.5) ≈ (0.4625 + 0.9 * 0.85) / 0.95 rtol = 1.0e-14
        @test complete_life_expectancy(qs, 0.5, UniformDeaths()) === complete_life_expectancy(qs, 0.5)
        for dd in (UniformDeaths(), ConstantForce(), Balducci()), age in (0.25, 0.5, 1.75)
            # the integral of survival conditional on surviving to `age`, as `survival` gives it
            ref = quadgk(u -> survival(qs, age, u, dd), age, ceil(age), 2)[1]
            @test complete_life_expectancy(qs, age, dd) ≈ ref rtol = 1.0e-12
            # survival from a fractional age multiplies with the expectancy at the next whole age
            s = survival(qs, age, ceil(age), dd)
            @test complete_life_expectancy(qs, age, dd) ≈
                quadgk(u -> survival(qs, age, u, dd), age, ceil(age))[1] + s * complete_life_expectancy(qs, ceil(Int, age), dd) rtol = 1.0e-14
        end
        # in the last year before the last age, only the part to the last age counts
        @test complete_life_expectancy(qs, 1.5) ≈ quadgk(u -> survival(qs, 1.5, u, UniformDeaths()), 1.5, 2)[1] rtol = 1.0e-14
        # Float32 rates and ages stay Float32
        q32 = UltimateMortality(Float32[0.1, 0.3, 1.0])
        @test (@inferred complete_life_expectancy(q32, 0.5f0)) isa Float32
        @test (@inferred complete_life_expectancy(qs, 0.5, ConstantForce())) isa Float64
        @test (@inferred complete_life_expectancy(qs, 1)) isa Float64
        @test complete_life_expectancy(q32, 0.5f0) ≈ complete_life_expectancy(qs, 0.5) rtol = 1.0e-6
    end

end
