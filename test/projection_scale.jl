@testset "projection scales" begin
    # inspired by Statutory Valuation of Life Insurance Liabilities, vol 5 chapter 21
    @testset "basic values" begin
        imp = MortalityTables.table("Projection Scale G2 – Male, ANB")
        @test imp[0] == 0.01
        @test imp[105] == 0.000
        @test imp[65] == 0.015
    end

    @testset "improved mortality" begin
        imp = MortalityTables.table("Projection Scale G2 – Male, ANB")
        mort = MortalityTables.table("2012 IAM Period Table – Male, ANB").ultimate

        @test mort[66] * prod(1 .- imp[65:65]) ≈ 8.41978 / 1000
        @test isapprox(mort[85] * prod(1 .- imp[65:84]), 44.60138 / 1000, atol=1e-5)

    end
end 