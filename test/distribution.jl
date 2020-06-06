# Reference: Experience Study Calculations, 2016, Society of Actuaries
# https://www.soa.org/globalassets/assets/Files/Research/2016-10-experience-study-calculations.pdf
# examples on pages 41-45
@testset "Fractional Year and distribution of deaths" begin
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
            @test survivorship(soa_mort,0,1,method) == 0.88
            @test survivorship(soa_mort,1,method) == 0.88
            @test cumulative_decrement(soa_mort,0,1,method) == 0.12
            @test cumulative_decrement(soa_mort,1,method) == 0.12

            # test floating point whole ages
            @test cumulative_decrement(soa_mort,0,12/12,method) == 0.12
            @test cumulative_decrement(soa_mort,12/12,method) == 0.12
        end
    
        # test fractional ages
        for i = 1:length(methods)
            for (t, target) in time_targets
                @test round(survivorship(soa_mort,t, methods[i]), digits = 4) ==
                      target[i]
                @test round(cumulative_decrement(soa_mort,t, methods[i]), digits = 4) ==
                      round(1 - target[i], digits = 4)
            end
        end

        # test time zero when given distribution of deaths
        for m in methods
            @test cumulative_decrement(soa_mort,0,m) == 0.0
            @test survivorship(soa_mort,0,m) == 1.0
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
            @test cumulative_decrement(mort,1,method) ≈ 0.20
            @test cumulative_decrement(mort,2,method) ≈ 0.20 + 0.8 * 0.5
            @test cumulative_decrement(mort,1,2,method) ≈ 0.50
        end
    
        # test fractional ages
        for i = 1:length(methods)
            for (t, target) in time_targets
                @test round(survivorship(mort, t, methods[i]), digits = 4) ==
                      target[i]
                @test round(cumulative_decrement(mort,t, methods[i]), digits = 4) ==
                      round(1 - target[i], digits = 4)
            end
        end
    end

    @testset "Error when asking for a fractional without assumption" begin
        mort = UltimateMortality([0.20, 0.50])
        @test_throws MethodError cumulative_decrement(mort, 0, 1, 1.5)
        @test_throws MethodError survivorship(mort, 0, 1, 1.5)
    end
end

