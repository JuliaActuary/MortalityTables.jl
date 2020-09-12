@testset "Parameterized Models" begin

    @testset "Makeham" begin

        g = Gompertz(a=0.0002,b=.13)

        @test survivorship(g,45) ≈ 0.5870365016720939
        @test survivorship(g,45,46) ≈ 0.9285202788707242
        @test decrement(g,45,46) ≈ 1 - 0.9285202788707242
        @test hazard(g,45) ≈ 0.06944687609574696

        @testset "AMLCR" begin
            m = Makeham(a=2.7e-6,b=log(1.124), c= 0.00022)

            @test MortalityTables.μ(m, 20) == 0.00022 + 2.7e-6 * 1.124^20
            @test m[20] == MortalityTables.μ(m, 20)
            @test m(20) == MortalityTables.μ(m, 20)
            
            # vs manually calculated (via QuadGK) integrals
            @test decrement(m, 20, 25) ≈ 0.0012891622754368504
            @test survivorship(m, 20, 25) ≈ 1 - 0.0012891622754368504
            @test decrement(m, 25) ≈ 0.005888764668801838
            @test survivorship(m, 25) ≈ 1 - 0.005888764668801838

            # these values come from the 'Standard Select and Ultimate Survival Model'
            # from Actuarial Mathematics for Life Contingent Risks, 2nd end

            ℓ = 100_000
            ℓs = [survivorship(m, 20, age) for age in 21:100] .* ℓ

            ℓ_age(x) = ℓs[x - 20]
            @test isapprox(ℓ_age(21),  99975.04, atol = 0.01)
            @test isapprox(ℓ_age(31),  99695.83, atol = 0.01)
            @test isapprox(ℓ_age(82),  70507.19, atol = 0.01)
            @test isapprox(ℓ_age(100),  6248.17, atol = 0.01)
        end

    end

    @testset "Gompertz and Makeham equality" begin

        # Gompertz is Makeham's where c = 0
        m = Makeham( a=2.7e-6,b= 1.124,c=0.0)
        g = Gompertz(a=2.7e-6,b= 1.124)

        for age ∈ 20:100
            @test survivorship(m, age) == survivorship(g, age)
            @test survivorship(m, age, 1) == survivorship(g, age, 1)
        end
    end

end