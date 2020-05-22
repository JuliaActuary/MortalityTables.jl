@testset "Makeham" begin
    # these values come from the 'Standard 
    # from Actuarial Mathematics for Life Contingent Risks, 2nd end


    m = Makeham(0.00022,2.7e-6,1.124)

    ℓ = 100_000
    ℓs = cumprod([p(m,age) for age ∈ 20:100]) .* ℓ

    ℓ_age(x) = ℓs[x - 20]
    @test isapprox( ℓ_age(21) ,  99975.04, atol=0.01)
    @test isapprox( ℓ_age(31) ,  99695.83, atol=0.01)
    @test isapprox( ℓ_age(82) ,  70507.19, atol=0.01)
    @test isapprox( ℓ_age(100) ,  6248.17, atol=0.01)

end