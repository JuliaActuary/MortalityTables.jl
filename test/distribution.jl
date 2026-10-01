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

    @testset "fractional starting ages condition on survival to the start" begin
        # survival from age 0 to y, built from whole-year products and the
        # within-year fraction starting at a birthday (independent of `survival`)
        within(::Uniform, q, u) = 1 - u * q
        within(::Constant, q, u) = (1 - q)^u
        within(::Balducci, q, u) = (1 - q) / (1 - (1 - u) * q)
        function S0(v, y, dd)
            x = floor(Int, y)
            return prod(1 - v[k] for k in firstindex(v):(x - 1); init = 1.0) * (y == x ? 1.0 : within(dd, v[x], y - x))
        end

        q = UltimateMortality([0.2, 0.3, 1.0])
        @test survival(q, 0.5, 1, Uniform()) ≈ 0.8 / 0.9
        @test survival(q, 0.5, 1, Balducci()) ≈ 0.9
        @test survival(q, 0.5, 1, Constant()) ≈ sqrt(0.8)
        for dd in (Uniform(), Balducci(), Constant())
            @test survival(q, 0, 0.5, dd) * survival(q, 0.5, 1, dd) ≈ 0.8
        end

        # a select row whose select period (ages 40:42) runs into the ultimate rates
        ult = UltimateMortality([0.01 * k for k in 1:20], start_age = 40)
        row = MortalityTables._select_row(40, [0.005, 0.02, 0.07], ult)
        ages = [0.0, 0.25, 0.5, 1.0, 1.4, 2.75]
        for (v, offset) in ((q, 0), (row, 40)), dd in (Uniform(), Balducci(), Constant())
            points = offset .+ ages
            for (i, a) in enumerate(points), c in points[(i + 1):end]
                @test survival(v, a, c, dd) ≈ S0(v, c, dd) / S0(v, a, dd)
                @test decrement(v, a, c, dd) ≈ 1 - S0(v, c, dd) / S0(v, a, dd)
                for b in points
                    a < b < c || continue
                    @test survival(v, a, c, dd) ≈ survival(v, a, b, dd) * survival(v, b, c, dd)
                end
            end
        end
        # crossing from the select period (age 42.5) into the ultimate rates (age 43.5)
        @test survival(row, 42.5, 43.5, Uniform()) ≈ (1 - 0.07) / (1 - 0.5 * 0.07) * (1 - 0.5 * ult[43])
    end
 
end

