@testset "Parameterized Models" begin
    @testset "Makeham" begin


        m = Makeham(0.00022,2.7e-6,1.124)

        @test MortalityTables.μ(m,20) == 0.00022 + 2.7e-6 * 1.124 ^ 20
        @test m[20] == MortalityTables.μ(m,20)
        @test m(20) == MortalityTables.μ(m,20)
        
        # vs manually calculated (via QuadGK) integrals
        @test cumulative_decrement(m,20,25) ≈ 0.0012891622754368504
        @test survivorship(m,20,25) ≈ 1 - 0.0012891622754368504
        @test cumulative_decrement(m,25) ≈ 0.005888764668801838
        @test survivorship(m,25) ≈ 1 - 0.005888764668801838

        # these values come from the 'Standard Select and Ultimate Survival Model'
        # from Actuarial Mathematics for Life Contingent Risks, 2nd end

        ℓ = 100_000
        ℓs = [survivorship(m,20,age) for age in 21:100] .* ℓ

        ℓ_age(x) = ℓs[x - 20]
        @test isapprox( ℓ_age(21) ,  99975.04, atol=0.01)
        @test isapprox( ℓ_age(31) ,  99695.83, atol=0.01)
        @test isapprox( ℓ_age(82) ,  70507.19, atol=0.01)
        @test isapprox( ℓ_age(100) ,  6248.17, atol=0.01)

    end

    @testset "Gompertz" begin

        # Gompertz is Makeham's where a = 0
        m = Makeham(0.0,2.7e-6,1.124)
        g = Gompertz(2.7e-6,1.124)

        for age ∈ 20:100
            @test survivorship(m,age) == survivorship(g,age)
            @test survivorship(m,age,1) == survivorship(g,age,1)
        end
    end
end