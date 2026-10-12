@testset "age basis conversion" begin
    @testset "published 1980 CSO values" begin
        nearest = MortalityTables.table(42).ultimate
        published_last = MortalityTables.table(41).ultimate
        converted = age_nearest_to_age_last(nearest)

        @test axes(converted) == axes(nearest) == axes(published_last)
        @test all(isapprox.(converted, published_last; atol=5.1e-6))
    end

    @testset "inverse and age indexes" begin
        nearest = UltimateMortality(Float32[0.01, 0.02, 0.04, 0.10, 1.0]; start_age = 40)
        last_birthday = age_nearest_to_age_last(nearest)
        round_trip = age_last_to_age_nearest(last_birthday)

        @test axes(last_birthday) == (40:44,)
        @test eltype(last_birthday) == Float32
        @test round_trip ≈ nearest
    end

    @testset "round trip through a zero rate" begin
        for T in (Float32, Float64)
            rates = T[0.0, 0.2, 1.0]
            round_trip = age_last_to_age_nearest(age_nearest_to_age_last(rates))
            @test round_trip ≈ rates
            @test iszero(first(round_trip))
            @test_throws ArgumentError age_last_to_age_nearest(T[0.001, 0.6, 1.0])
        end
    end

    @testset "invalid inputs" begin
        @test_throws ArgumentError age_nearest_to_age_last(Float64[])
        @test_throws ArgumentError age_nearest_to_age_last([0.1, 0.2])
        @test_throws ArgumentError age_nearest_to_age_last([0.1, 1.1, 1.0])
        @test_throws ArgumentError age_last_to_age_nearest([0.0, 1.0])
    end
end
