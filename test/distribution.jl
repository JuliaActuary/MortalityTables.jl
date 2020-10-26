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
            @test survival(soa_mort, 0, 1, method) == 0.88
            @test survival(soa_mort, 1, method) == 0.88
            @test decrement(soa_mort, 0, 1, method) == 0.12
            @test decrement(soa_mort, 1, method) == 0.12

            # test floating point whole ages
            @test decrement(soa_mort, 0, 12 / 12, method) == 0.12
            @test decrement(soa_mort, 12 / 12, method) == 0.12
        end
    
        # test fractional ages
        for i = 1:length(methods)
            for (t, target) in time_targets
                @test round(survival(soa_mort, t, methods[i]), digits=4) ==
                      target[i]

                @test round(MortalityTables.decrement_partial_year(soa_mort, 0, t, methods[i]),
                    digits=4) ≈  1 - target[i]

                @test round(decrement(soa_mort, t, methods[i]), digits=4) ==
                      round(1 - target[i], digits=4)
            end
        end

        # test time zero when given distribution of deaths
        for m in methods
            @test decrement(soa_mort, 0, m) == 0.0
            @test survival(soa_mort, 0, m) == 1.0
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
            @test decrement(mort, 1, method) ≈ 0.20
            @test decrement(mort, 2, method) ≈ 0.20 + 0.8 * 0.5
            @test decrement(mort, 1, 2, method) ≈ 0.50
        end
    
        # test fractional ages
        for i = 1:length(methods)
            for (t, target) in time_targets
                @test round(survival(mort, t, methods[i]), digits=4) ==
                      target[i]
                @test round(decrement(mort, t, methods[i]), digits=4) ==
                      round(1 - target[i], digits=4)
            end
        end
    end

    @testset "Error when asking for a fractional without assumption" begin
        mort = UltimateMortality([0.20, 0.50])
        @test_throws MethodError decrement(mort, 0, 1, 1.5)
        @test_throws MethodError survival(mort, 0, 1, 1.5)
    end

    @testset "Issue #60 - starting with non-integer ages" begin
        m = UltimateMortality([0.5 for i in 1:8])
        
        @test survival(m, 1, 2) ≈ 0.5
        @test survival(m, 1.5, 2.5, Constant()) ≈ 0.5
        @test survival(m, 1.5, 3.5, Constant()) ≈ 0.25
        @test survival(m, 1.5, 1.5 + eps(), Constant()) ≈ 1.0
        @test survival(m, 1, 1 + eps(), Constant()) ≈ 1.0

    end

    @testset "Issue #88 - Stackoverflow" begin
        @test survival(MortalityTables.mortality_vector([0.5,0.5],start_age=50),50,50.5,MortalityTables.Uniform()) ≈ 1 - 0.5 * 0.5
    end
 
end

