
@testset "basic MortalityTable" begin
    @testset "basic structure" begin
        q1 = UltimateMortality([i for i = 0:19], start_age = 0)
        @test q1[0] == 0
        @test q1[1] == 1
        @test q1[0:1] == [0, 1]
        @test omega(q1) == 19
        @test_throws BoundsError q1[omega(q1)+1]

        # non-zero start age
        q2 = UltimateMortality([i for i = 0:9], start_age = 5)
        @test_throws BoundsError q2[4]
        @test q2[5] == 0

        # select strucutre
        s = [ia + d for ia = 0:5, d = 0:4]

        q3 = SelectMortality(s, q1, start_age = 0)
        @test q3[0][0] == 0
        @test q3[0][0:1] == [0, 1]
        @test omega(q3[0]) == 19
        @test_throws BoundsError q3[omega(q3[0])+1]
        @test q3[0][19] == 19

        @test q3[5][5] == 5


        mt1 = MortalityTable(q3, q1)

        @test mt1.select[0][1] == 1
        @test mt1.ultimate[1] == 1

        mt2 = MortalityTable(q1)

        @test mt2.ultimate[0] == 0
        @test_throws MethodError mt2[0]  # tables are opaque containers

    end

    @testset "_select_row" begin
        ult = UltimateMortality([i / 100 for i = 0:100])

        # (a) leading missing: the row is still indexed from the issue age
        row = MortalityTables._select_row(15, [missing, 0.2, 0.3], ult)
        @test firstindex(row) == 15
        @test row[15] === missing
        @test row[16] == 0.2
        @test row[17] == 0.3
        @test row[18] == ult[18]
        @test omega(row) == 100
        @test eltype(row) == Union{Missing,Float64}

        # a row without missing keeps a Float64 element type
        row = MortalityTables._select_row(15, [0.1, 0.2, 0.3], ult)
        @test eltype(row) == Float64
        @test row[15:17] == [0.1, 0.2, 0.3]
        @test row[18:100] == ult[18:100]

        # (b) select period ending exactly at omega: no ultimate tail
        row = MortalityTables._select_row(95, fill(0.5, 6), ult)
        @test omega(row) == 100
        @test all(row .== 0.5)

        # (c) select period extending past omega: the select rates alone
        row = MortalityTables._select_row(95, fill(0.5, 10), ult)
        @test omega(row) == 104
        @test length(row) == 10
        @test all(row .== 0.5)
    end

    @testset "off-aligned select and ult" begin
        select_matrix = [(i + j - 1) / 100 for i = 0:10, j = 1:20]
        ult = UltimateMortality([i / 100 for i = 18:100], start_age = 18)

        select = SelectMortality(select_matrix, ult)

        for issue_age = 0:10
            for dur = 1:30
                @test select[issue_age][issue_age+dur-1] == (issue_age + dur - 1) / 100
            end
        end
    end

    # test time zero accumlated force
    @testset "accumulated force" begin
        q4 = UltimateMortality([0.1, 0.3, 0.6, 1])

        @test survival(q4, 0) ≈ 1
        @test decrement(q4, 0) ≈ 0

        @test survival(q4, 1) ≈ 0.9
        @test decrement(q4, 1) ≈ 0.1

        @test survival(q4, 1, 1) ≈ 1.0
        @test survival(q4, 1, 2) ≈ 0.7
        @test decrement(q4, 1, 1) ≈ 0.0
        @test decrement(q4, 1, 2) ≈ 0.3

        @test survival(q4, 1, 4) ≈ 0.0
        @test decrement(q4, 1, 4) ≈ 1.0


        # the table has no rates before age 0
        @test_throws BoundsError survival(q4, -1)
        @test_throws BoundsError survival(q4, 4, -1)
        @test_throws BoundsError decrement(q4, -1)
        @test_throws BoundsError decrement(q4, 4, -1)
    end

    @testset "small decrements keep their precision" begin
        q = UltimateMortality([1e-18, 1e-6, 0.3, 0.0, 1.0])
        # 1 - (1 - q) rounds 1e-18 to zero and 1e-6 to about 11 correct digits
        @test decrement(q, 0, 1) == 1e-18
        @test decrement(q, 1, 2) == 1e-6
        @test decrement(q, 0, 2) == 1.0000000000009999e-6   # 1e-18 + 1e-6 - 1e-24
        @test decrement(q, 0, 3) ≈ 1 - survival(q, 0, 3)
        # against a high-precision table, whole and fractional ages
        qbig = UltimateMortality(big.([1e-18, 1e-6, 0.3, 0.0, 1.0]))
        for dd in (Uniform(), Balducci(), Constant()), (a, b) in ((0.25, 0.75), (0.5, 1.5), (0.0, 1.25), (1.1, 1.9), (0.3, 2.6))
            @test decrement(q, a, b, dd) ≈ Float64(decrement(qbig, big(a), big(b), dd)) rtol = 1e-13
            @test decrement(q, a, b, dd) ≈ 1 - survival(q, a, b, dd) atol = 1e-15
        end
        @test decrement(q, 0.25, 0.75, Constant()) > 0   # not rounded away
        @test decrement(q, 1, 2, Uniform()) == 1e-6

        # terminal rates: q = 0 contributes nothing, q = 1 exhausts survival
        @test decrement(q, 3, 4) == 0.0
        @test survival(q, 3, 4) == 1.0
        @test decrement(q, 4, 5) == 1.0
        @test survival(q, 4, 5) == 0.0
        @test decrement(q, 2, 5) == 1.0
        for dd in (Uniform(), Balducci(), Constant())
            # equal endpoints are exactly the identity, even at a rate of one
            for a in (2.0, 3.5, 4.0, 4.5)
                @test survival(q, a, a, dd) == 1.0
                @test decrement(q, a, a, dd) == 0.0
            end
            # a piece starting at a birthday, and one ending at the next birthday
            @test decrement(q, 2, 2.5, dd) ≈ 1 - survival(q, 2, 2.5, dd)
            @test decrement(q, 2.5, 3, dd) ≈ 1 - survival(q, 2.5, 3, dd)
            @test decrement(q, 4.0, 4.5, dd) == (dd isa Uniform ? 0.5 : 1.0)
            @test decrement(q, 3.25, 3.75, dd) == 0.0
        end
        # Uniform and Balducci agree with q at s = 0, t = 1 of the year
        for dd in (Uniform(), Balducci())
            @test MortalityTables.decrement_partial_year(q, 2, 3, dd) == 0.3
            @test decrement(q, 2.0, 3.0, dd) == 0.3
        end
        # under Uniform, a rate of one leaves positive survival inside the year
        @test survival(q, 4.0, 4.5, Uniform()) == 0.5
        # beyond a rate of one under Constant or Balducci there is no survival to condition on:
        # the formulas are evaluated as written and stay finite
        for dd in (Balducci(), Constant())
            @test survival(q, 4.0, 4.3, dd) == 0.0
            @test isfinite(survival(q, 4.3, 4.6, dd))
            @test isfinite(decrement(q, 4.3, 4.6, dd))
        end
        @test survival(q, 4.3, 4.6, Balducci()) ≈ 0.5   # (1 - (1 - s)q) / (1 - (1 - t)q) at q = 1
    end

    @testset "results take the rates' numeric type, empty intervals included" begin
        q64 = UltimateMortality([0.1, 0.3, 1.0])
        q32 = UltimateMortality(Float32[0.1, 0.3, 1.0])
        qbig = UltimateMortality(big.([0.1, 0.3, 1.0]))
        for (q, T) in ((q64, Float64), (q32, Float32), (qbig, BigFloat))
            # whole ages: the rates' type
            @test survival(q, 1, 1) isa T
            @test survival(q, 1, 1) == 1
            @test survival(q, 0, 2) isa T
            @test survival(q, 2, 0) isa T
            @test decrement(q, 1, 1) isa T
            @test decrement(q, 1, 1) == 0
            @test decrement(q, 0, 2) isa T
            @test decrement(q, 2, 0) isa T
            @test life_expectancy(q, 2) isa T
            @test life_expectancy(q, 2) == 0
            @test life_expectancy(q, 0) isa T
            @test life_expectancy(q, 0) ≈ 0.9 + 0.9 * 0.7
            # fractional ages promote with the ages' type
            S = promote_type(T, Float64)
            for dd in (Uniform(), Balducci(), Constant())
                @test survival(q, 0.5, 0.5, dd) isa S
                @test survival(q, 0.5, 1.5, dd) isa S
                @test decrement(q, 0.5, 0.5, dd) isa S
                @test decrement(q, 0.5, 1.5, dd) isa S
            end
            @test life_expectancy(q, 2, Uniform()) isa S
            @test life_expectancy(q, 2, Uniform()) == 0
            @test life_expectancy(q, 0, Uniform()) isa S
        end
        # Float32 rates and Float32 ages stay Float32
        @test survival(q32, 0.5f0, 1.5f0, Uniform()) isa Float32
        @test decrement(q32, 0.5f0, 1.5f0, Constant()) isa Float32
        # rates that admit `missing` (a select row with a gap) are typed by their numbers
        row = MortalityTables._select_row(0, [missing, 0.2], q64)
        @test survival(row, 1, 1) isa Float64
        @test survival(row, 1, 1) == 1
        @test survival(row, 1, 2) ≈ 0.8
        @test decrement(row, 1, 1) === 0.0
        # a missing rate at the last age is never used: the expectancy there is zero, and one
        # year earlier it is the survival to the last age
        gap = UltimateMortality([0.1, missing])
        @test life_expectancy(gap, 1) === 0.0
        @test life_expectancy(gap, 0) ≈ 0.9
        # a missing rate that is used propagates
        @test ismissing(life_expectancy(UltimateMortality([missing, 0.2, 0.3]), 0))
    end

    @testset "reversed intervals are reverse factors" begin
        q4 = UltimateMortality([0.1, 0.3, 0.6, 1])
        # the inverse of the forward survival, and its decrement 1 - 1/S is negative
        @test survival(q4, 2, 0) ≈ 1 / (0.9 * 0.7)
        @test decrement(q4, 2, 0) ≈ 1 - 1 / (0.9 * 0.7)
        @test decrement(q4, 2, 0) < 0
        @test survival(q4, 1, 0) * survival(q4, 0, 1) ≈ 1
        # a zero forward survival has no finite inverse
        @test survival(q4, 4, 0) == Inf
        for dd in (Uniform(), Balducci(), Constant())
            @test survival(q4, 2.5, 0.25, dd) ≈ 1 / survival(q4, 0.25, 2.5, dd)
            @test survival(q4, 1.75, 1.25, dd) ≈ 1 / survival(q4, 1.25, 1.75, dd)   # within one year
            @test decrement(q4, 2.5, 0.25, dd) ≈ 1 - 1 / survival(q4, 0.25, 2.5, dd)
        end

        # survival composes over any three ages, whatever their order
        ult = UltimateMortality([0.01 * k for k in 1:20], start_age = 40)
        row = MortalityTables._select_row(40, [0.005, 0.02, 0.07], ult)
        for v in (ult, row), dd in (Uniform(), Balducci(), Constant())
            ages = (40.0, 40.3, 41.0, 42.5, 43.75, 45.0)
            for a in ages, b in ages, c in ages
                @test survival(v, a, c, dd) ≈ survival(v, a, b, dd) * survival(v, b, c, dd)
            end
            for a in 40:45, b in 40:45
                @test survival(v, a, b) * survival(v, b, a) ≈ 1
            end
        end

        # projecting a population backward: 1,000 lives at 45 imply the number expected at 40
        l45 = 1000.0
        l40 = l45 * survival(ult, 45, 40)
        @test l40 ≈ l45 / prod(1 - ult[x] for x in 40:44)
        @test l40 * survival(ult, 40, 45) ≈ l45
    end

    @testset "Metadata" begin
        d = TableMetaData()

        @test isnothing(d.name)

        d = TableMetaData(name = "test")
        @test d.name == "test"
    end

    @testset "show" begin
        # a custom table has no mort.SOA.org id, so no link should be printed
        s = sprint(show, MIME"text/plain"(), MortalityTable(UltimateMortality([0.1, 0.2])))
        @test !occursin("TableIdentity=nothing", s)
        @test !occursin("mort.SOA.org ID", s)
        @test occursin("Description", s)

        s = sprint(show, MIME"text/plain"(), MortalityTable(UltimateMortality([0.1, 0.2]), metadata = TableMetaData(id = "42")))
        @test occursin("TableIdentity=42", s)
    end

    @testset "mortality_vector" begin
        v = [i for i = 3:10]
        q = mortality_vector(v, start_age = 3)
        @test q[3] == 3
        @test q[10] == 10

        q = mortality_vector(collect(0:5))
        @test q[0] == 0
        @test q[5] == 5

        # mortality_vector is an alias of UltimateMortality
        @test mortality_vector(v, start_age = 3) == UltimateMortality(v, start_age = 3)
    end

    @testset "UltimateMortality accepts any AbstractVector" begin
        # a range
        r = UltimateMortality(0:0.1:1)
        @test r[0] == 0.0
        @test r[10] == 1.0

        # a view
        base = [0.1, 0.2, 0.3, 0.4]
        vw = UltimateMortality(view(base, 2:4), start_age = 1)
        @test vw[1] == 0.2
        @test vw[3] == 0.4

        # a vector containing missing (as the XTbML and CSV loaders produce)
        mv = UltimateMortality([0.1, missing])
        @test mv[0] == 0.1
        @test mv[1] === missing

        # no copy is made
        q = UltimateMortality(base)
        base[1] = 0.5
        @test q[0] == 0.5

        # start_age sets the first age whatever the input's own axes: a vector already indexed
        # by age is re-anchored, not shifted relative to its old first age
        for v in ([0.1, 0.2], view([0.0, 0.1, 0.2], 2:3), UltimateMortality([0.1, 0.2]; start_age = 40),
                  UltimateMortality([0.1, 0.2]; start_age = -5))
            m = UltimateMortality(v; start_age = 7)
            @test axes(m, 1) == 7:8
            @test m[7] == 0.1 && m[8] == 0.2
        end
    end

    @testset "utility functions" begin
        @test MortalityTables._decrement(1.0, 0.05) ≈ 1.0 * (1 - 0.05)
    end

end
