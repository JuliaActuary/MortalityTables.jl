using FinanceCore: AbstractDeflator, factor, intensity, compose, Continuous, Cashflow, pv
using MortalityTables: ParametricMortality, AbstractMortality

@testset "AbstractDeflator integration" begin
    @testset "Type hierarchy" begin
        @test AbstractMortality <: AbstractDeflator
        @test ParametricMortality <: AbstractMortality
        @test MortalityVector <: AbstractMortality
        @test ShiftedMortality <: AbstractMortality
    end

    @testset "UltimateMortality returns MortalityVector" begin
        rates = [x / 100 for x in 0:99]
        m = UltimateMortality(rates)
        @test m isa MortalityVector
        @test m isa AbstractDeflator
        # AbstractArray-like interface still works (back-compat)
        @test m[5] == 0.05
        @test length(m) == 100
        @test firstindex(m) == 0
        @test lastindex(m) == 99
        @test parent(m) isa AbstractArray  # underlying OffsetArray

        # start_age option
        m18 = UltimateMortality(rates, start_age = 18)
        @test firstindex(m18) == 18
        @test m18[18] == 0.0
    end

    @testset "SelectMortality returns MortalityVectors per issue age" begin
        ult = UltimateMortality([x / 100 for x in 0:99])
        select_matrix = rand(50, 10)
        sel = SelectMortality(select_matrix, ult; start_age = 0)
        # Indexing by issue age returns a MortalityVector
        @test sel[0] isa MortalityVector
        # Underlying array access still works
        @test sel[0][0] == select_matrix[1, 1]
    end

    @testset "factor/intensity on ParametricMortality" begin
        m = Makeham(a = 0.0002, b = 0.13, c = 0.001)
        @test m isa AbstractMortality
        @test m isa AbstractDeflator
        @test factor(m, 65) ≈ survival(m, 65)
        @test factor(m, 65, 70) ≈ survival(m, 65, 70)
        @test intensity(m, 65) == hazard(m, 65)
    end

    @testset "factor on MortalityVector" begin
        rates = [x / 100 for x in 0:99]
        m = UltimateMortality(rates)
        # Two-arg: integer ages
        @test factor(m, 65, 70) ≈ survival(m, 65, 70)
        @test factor(m, 65, 70) ≈ prod(1 - rates[k+1] for k in 65:69)
        # Single-arg: anchored at firstindex (= 0 here, so survival from age 0)
        @test factor(m, 5) ≈ survival(m, 0, 5)
        # Fractional time without DD raises
        @test_throws ArgumentError factor(m, 5.5)
        @test_throws ArgumentError factor(m, 0.5, 5.0)
        # Fractional time with DD works
        @test factor(m, 65.0, 70.0, Uniform()) ≈ survival(m, 65.0, 70.0, Uniform())
    end

    @testset "at_age: align to valuation age" begin
        # Parametric model
        m = Makeham(a = 2.5e-5, b = 0.10, c = 1e-4)
        m_at_65 = at_age(m, 65)
        @test m_at_65 isa ShiftedMortality
        @test m_at_65 isa AbstractMortality

        # factor on years-from-valuation axis
        @test factor(m_at_65, 5) ≈ survival(m, 65, 70)
        @test factor(m_at_65, 0, 5) ≈ survival(m, 65, 70)
        @test intensity(m_at_65, 0) ≈ hazard(m, 65)
        @test intensity(m_at_65, 5) ≈ hazard(m, 70)

        # MortalityVector — at_age shifts the origin
        rates = [x / 100 for x in 0:99]
        mv = UltimateMortality(rates)
        mv_at_65 = at_age(mv, 65)
        @test factor(mv_at_65, 5) ≈ factor(mv, 65, 70)
        @test factor(mv_at_65, 0, 5) ≈ factor(mv, 65, 70)
    end

    @testset "compose: yield × at_age(parametric)" begin
        yield = Continuous(0.03)
        m = Makeham(a = 2.5e-5, b = 0.10, c = 1e-4)
        deflator = compose(yield, at_age(m, 65))

        # 5-year EPV of $1 contingent on a 65-year-old surviving
        expected = exp(-0.03 * 5) * survival(m, 65, 70)
        @test pv(deflator, [Cashflow(1.0, 5)]) ≈ expected
    end

    @testset "compose: yield × at_age(MortalityVector)" begin
        yield = Continuous(0.03)
        rates = [x / 100 for x in 0:99]
        mv = UltimateMortality(rates)
        mort_at_65 = at_age(mv, 65)
        deflator = compose(yield, mort_at_65)

        # PV of $1 paid in 5 years contingent on the 65-year-old surviving
        expected = exp(-0.03 * 5) * survival(mv, 65, 70)
        @test factor(deflator, 0, 5) ≈ expected
        @test pv(deflator, [Cashflow(1.0, 5)]) ≈ expected
    end

    @testset "MortalityTable container is NOT a deflator (extract a vector first)" begin
        # The MortalityTable abstract type stays a metadata container.
        @test !(MortalityTable <: AbstractDeflator)
    end

    @testset "survival/decrement aliases on AbstractMortality" begin
        m = Makeham(a = 2.5e-5, b = 0.10, c = 1e-4)
        @test survival(m, 65, 70) ≈ factor(m, 65, 70)
        @test decrement(m, 65, 70) ≈ 1 - factor(m, 65, 70)

        rates = [x / 100 for x in 0:99]
        mv = UltimateMortality(rates)
        @test survival(mv, 65, 70) ≈ factor(mv, 65, 70)
        @test decrement(mv, 65, 70) ≈ 1 - factor(mv, 65, 70)
    end
end
