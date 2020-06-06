@testset "Parameterized Models" begin
    @testset "Makeham" begin
        # these values come from the 'Standard 
        # from Actuarial Mathematics for Life Contingent Risks, 2nd end


        m = Makeham(0.00022,2.7e-6,1.124)

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