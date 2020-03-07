using MortalityTables
using Test

include("basic.jl")
include("XTbML.jl")

#load tables to be used in subsequent tests
tables = MortalityTables.tables()
@test length(tables) > 0


@testset "Mortality Assumption" begin

    vbt2001 = tables["2001 VBT Residual Standard Select and Ultimate - Male Nonsmoker, ANB"]


    # q(b,20,0.5) ≈ 0.00077 * # Balducci

    # uniform
    u = MortalityAssumption(vbt2001,Uniform())
    @test q(u,25,1) ≈ 0.00101
    @test p(u,25,0.5) ≈ 1 - 0.00101 * 0.5
    @test p(u,25,1.5) ≈ (1 - 0.00101) * (1 - 0.00104 * 0.5)

    #constant
    c = MortalityAssumption(vbt2001,Constant())
    @test q(c,25,1) ≈ 0.00101
    @test p(c,25,0.5) ≈ (1 - 0.00101) ^ 0.5
    @test p(c,25,1.5) ≈ (1 - 0.00101) * (1 - 0.00104) ^ 0.5

end
include("MortalityTables.jl")
