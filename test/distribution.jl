# Reference: Experience Study Calculations, 2016, Society of Actuaries
# https://www.soa.org/globalassets/assets/Files/Research/2016-10-experience-study-calculations.pdf
# examples on pages 41-45
@testset "Fractional Year and distribution of deaths" begin

    @testset "ALMCR §9.4" begin
        ℓ =[43302,42854,42081,41351,40050]

        ps =  ℓ ./ ℓ[1]
        qs = [1 - ps[t] / ps[t-1] for t in 2:5 ] 
        
        m1 = UltimateMortality(qs,65)

        @test p(m1,65,1,1) == ℓ[2] / ℓ[1]
        @test q(m1,65,1,1) == 1 - ℓ[2] / ℓ[1]
        @test q(m1,65,1,1,Uniform()) == 1 - ℓ[2] / ℓ[1]

        @test p(m1,65,1,2) == ℓ[3] / ℓ[1]
        @test q(m1,65,1,2) == 1 - ℓ[3] / ℓ[1]
        @test q(m1,65,1,2,Uniform()) == 1 - ℓ[3] / ℓ[1]

    end

    @testset "SOA Experience Study Calcuations distribution examples" begin
        soa_mort = UltimateMortality([0.12])

        methods = [Balducci(), Uniform(), Constant()]
        time_targets = Dict(
            1 / 12 => [0.9888, 0.9900, 0.9894],
            6 / 12 => [0.9362, 0.9400, 0.9381],
            7 / 12 => [0.9263, 0.9300, 0.9281],
            12 / 12 => [0.8800, 0.8800, 0.8800],
        )

        # test whole ages with assumption argument
        for method in methods
            @test q(soa_mort,0,1,1,method) == 0.12
        end
    
        # test fractional ages
        for i = 1:length(methods)
            for (t, target) in time_targets
                @test round(p(soa_mort, 0, 1, t, methods[i]), digits = 4) ==
                      target[i]
                @test round(q(soa_mort, 0, 1, t, methods[i]), digits = 4) ==
                      round(1 - target[i], digits = 4)
            end
        end
    end
    @testset "Multi-year examples" begin
        mort = UltimateMortality([0.20, 0.50])
        methods = [Balducci(), Uniform(), Constant()]

        # these sample values calculated manually
        time_targets = Dict(
            1 + 1 / 12 => [0.7385, 0.7667, 0.7551],
            1 + 6 / 12 => [0.5333, 0.6000, 0.5657],
            1 + 7 / 12 => [0.5053, 0.5667, 0.5339],
            1 + 12 / 12 => [0.4000, 0.4000, 0.4000],
        )
        
        # test whole ages with assumption argument
        for method in methods
            @test q(mort,0,1,1,method) ≈ 0.20
            @test q(mort,1,1,1,method) ≈ 0.50
        end
    
        # test fractional ages
        for i = 1:length(methods)
            for (t, target) in time_targets
                @test round(p(mort, 0, 1, t, methods[i]), digits = 4) ==
                      target[i]
                @test round(q(mort, 0, 1, t, methods[i]), digits = 4) ==
                      round(1 - target[i], digits = 4)
            end
        end
    end

    @testset "Error when asking for a fractional without assumption" begin
        mort = UltimateMortality([0.20, 0.50])
        @test_throws ArgumentError q(mort, 0, 1, 1.5)
        @test_throws ArgumentError p(mort, 0, 1, 1.5)
    end
end
