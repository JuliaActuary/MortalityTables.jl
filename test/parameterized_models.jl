using InteractiveUtils: subtypes
using ForwardDiff

@testset "Parameterized Models" begin

    @testset "concrete field types" begin
        # the subtypes are UnionAlls (e.g. Makeham{T}), so iterate them directly
        laws = subtypes(MortalityTables.ParametricMortality)
        @test length(laws) == 25
        for L in laws
            m = L()
            @test isbits(m)
            @test m isa L{Float64}
            @test (@inferred hazard(m, 50.0)) isa Float64
        end
        # mixed keyword types are promoted
        @test Makeham(a=1, b=0.5, c=0) isa Makeham{Float64}
        @test Gompertz() isa Makeham{Float64}
    end

    @testset "cumhazard is the primitive for every law" begin
        for L in subtypes(MortalityTables.ParametricMortality)
            m = L()
            @test cumhazard(m, 50) isa Real
            @test survival(m, 0) ≈ 1
            # any closed-form cumhazard must agree with integrating the law's own hazard
            @test cumhazard(m, 50) ≈ quadgk(a -> hazard(m, a), 0, 50)[1] rtol = 1e-6
            # the ratio form is only defined when survival has not underflowed to zero
            # (the default Quadratic and VanderMaen parameters do); the difference of
            # cumhazards used by the package is well defined either way
            if survival(m, 40) > 0
                @test survival(m, 40, 60) ≈ survival(m, 60) / survival(m, 40)
            else
                @test 0 <= survival(m, 40, 60) <= 1
            end
            @test decrement(m, 40, 60) ≈ 1 - survival(m, 40, 60) atol = 1e-15
        end
    end

    @testset "zero growth: the constant-hazard limit" begin
        # b = 0 is a constant hazard: a + c for Makeham, a/(1 + a) for Kannisto. The closed
        # forms used to evaluate 0/0 there.
        for (m, λ) in ((Makeham(a = 0.001, b = 0.0, c = 0.002), 0.003), (Gompertz(a = 0.001, b = 0.0), 0.001),
                       (MortalityTables.Kannisto(a = 0.5, b = 0.0), 0.5 / 1.5))
            @test cumhazard(m, 50) ≈ 50λ rtol = 1e-14
            @test survival(m, 40, 41) ≈ exp(-λ) rtol = 1e-14
        end
        # continuous through b = 0, with the closed-form derivative in b there on both sides of
        # the series cutoff: ∂H/∂b = a·x²/2 (Makeham) and a·x²/(2(1 + a)²) (Kannisto)
        laws = ((b -> Makeham(a = 0.001, b = b, c = 0.002), 0.001 * 50^2 / 2),
                (b -> MortalityTables.Kannisto(a = 0.5, b = b), 0.5 * 50^2 / (2 * 1.5^2)))
        for (law, dH) in laws
            H(b) = cumhazard(law(b), 50)
            for b in (-1e-4, -1e-6, 1e-6, 1e-4)
                @test H(b) ≈ quadgk(x -> hazard(law(b), x), 0, 50)[1] rtol = 1e-10
            end
            for h in (1e-8, 1e-7, 1e-5)   # both sides of the |b·age| < 1e-5 series cutoff
                @test (H(h) - H(-h)) / 2h ≈ dH rtol = 1e-6
            end
        end
    end

    @testset "large growth and infinite ages" begin
        # A growing hazard exhausts survival; a decaying one has a finite total hazard.
        for m in (Makeham(), Gompertz(), MortalityTables.Kannisto())
            @test cumhazard(m, Inf) == Inf
            @test survival(m, Inf) == 0.0
        end
        # A zero coefficient contributes nothing, even at an infinite age: a constant hazard
        # (b = 0) still exhausts survival, and a law with no hazard keeps it.
        for m in (Makeham(a = 0.001, b = 0.0, c = 0.002), Gompertz(a = 0.001, b = 0.0),
                  MortalityTables.Kannisto(a = 0.5, b = 0.0), Makeham(a = 0.0, b = 0.13, c = 0.002))
            @test cumhazard(m, Inf) == Inf
            @test survival(m, Inf) == 0.0
        end
        for m in (Makeham(a = 0.0, b = 0.13, c = 0.0), Makeham(a = 0.0, b = 0.0, c = 0.0),
                  MortalityTables.Kannisto(a = 0.0, b = 0.13))
            @test cumhazard(m, Inf) == 0.0
            @test survival(m, Inf) == 1.0
        end
        @test cumhazard(Gompertz(a = 0.001, b = -0.1), Inf) ≈ 0.01 rtol = 1e-14
        @test cumhazard(MortalityTables.Kannisto(a = 0.5, b = -0.1), Inf) ≈ log1p(0.5) / 0.1 rtol = 1e-14
        # both sides of the switch to the factored forms at |b·age| = 1, against the closed
        # forms in high precision
        H_makeham(a, b, c, x) = a / b * expm1(b * x) + c * x
        H_kannisto(a, b, x) = log((1 + a * exp(b * x)) / (1 + a)) / b
        for b in (0.13, -0.13), age in (1 / 0.13 * (1 - 1e-6), 1 / 0.13 * (1 + 1e-6))
            @test cumhazard(Makeham(a = 0.0002, b = b, c = 0.001), age) ≈
                  H_makeham(big(0.0002), big(b), big(0.001), big(age)) rtol = 1e-14
            @test cumhazard(MortalityTables.Kannisto(a = 0.5, b = b), age) ≈
                  H_kannisto(big(0.5), big(b), big(age)) rtol = 1e-14
        end
        # a·exp(b·age) far from 1 + a: the factored form would round u to -1
        m = MortalityTables.Kannisto(a = 1e18, b = -0.5)
        ref = exp(H_kannisto(big(1e18), big(-0.5), big(80)) - H_kannisto(big(1e18), big(-0.5), big(100)))
        @test survival(m, 80, 100) ≈ ref rtol = 1e-12
    end

    @testset "ratio-form hazards stay finite where exp overflows" begin
        K = MortalityTables
        # the original ratio formulas, evaluated in high precision as references
        beard(a, b, k, x) = a * exp(b * x) / (1 + k * a * exp(b * x))
        gg(a, b, γ, x) = a * exp(b * x) / (1 + a * γ / b * (exp(b * x) - 1))
        mart(a, b, c, d, k, x) = (a * exp(b * x) + c) / (1 + d * exp(b * x)) + k * exp(b * x)
        refs = (
            (K.Beard(), x -> beard(big(0.002), big(0.13), big(1.0), x)),
            (K.Beard(k = 0.5, b = -0.1), x -> beard(big(0.002), big(-0.1), big(0.5), x)),
            (K.MakehamBeard(), x -> beard(big(0.002), big(0.13), big(1.0), x) + big(0.01)),
            (K.Kannisto(), x -> beard(big(0.5), big(0.13), big(1.0), x)),
            (K.KannistoMakeham(), x -> beard(big(0.5), big(0.13), big(1.0), x) + big(0.001)),
            (K.GammaGompertz(), x -> gg(big(0.002), big(0.13), big(1.0), x)),
            (K.GammaGompertz(b = -0.13), x -> gg(big(0.002), big(-0.13), big(1.0), x)),
            (K.Martinelle(), x -> mart(big(0.001), big(0.13), big(0.001), big(0.1), big(0.001), x)),
            (K.Martinelle(b = -0.13), x -> mart(big(0.001), big(-0.13), big(0.001), big(0.1), big(0.001), x)),
        )
        for (m, ref) in refs, age in (0.0, 0.5, 5.0, 7.7, 20.0, 50.0, 80.0, 110.0)
            @test hazard(m, age) ≈ ref(big(age)) rtol = 1e-14
        end
        # far past the point where exp(b·age) overflows, the hazard tends to its bound
        @test hazard(K.Beard(k = 0.5), 1e4) == 2.0
        @test hazard(K.MakehamBeard(), 1e4) ≈ 1.01
        @test hazard(K.Kannisto(), 1e4) == 1.0
        @test hazard(K.KannistoMakeham(), 1e4) ≈ 1.001
        @test hazard(K.GammaGompertz(), 1e4) ≈ 0.13
        @test hazard(K.Martinelle(k = 0.0), 1e4) ≈ 0.01
        @test hazard(K.Beard(b = -0.1), Inf) == 0.0
        @test hazard(K.Kannisto(a = 0.0), Inf) == 0.0
        # GammaGompertz at b = 0 is a / (1 + a·γ·age)
        @test hazard(K.GammaGompertz(b = 0.0), 50.0) ≈ 0.002 / (1 + 0.002 * 50) rtol = 1e-14
        # Kannisto's hazard is bounded by one, so survival over a year in the high-growth
        # tail is about exp(-1), not NaN
        H_kannisto(a, b, x) = log((1 + a * exp(b * x)) / (1 + a)) / b
        for (a, b) in ((0.5, 10.0), (1e308, 0.1))
            m = K.Kannisto(a = a, b = b)
            ref = exp(H_kannisto(big(a), big(b), big(80)) - H_kannisto(big(a), big(b), big(81)))
            @test survival(m, 80, 81) ≈ ref rtol = 1e-12
        end
        @test survival(K.Kannisto(a = 0.5, b = 10.0), 80, 81) ≈ 0.36787944117144233 rtol = 1e-12
        @test isfinite(cumhazard(K.Kannisto(a = 1e308, b = 0.1), 5.0))   # series branch
    end

    @testset "ratio-form hazards follow the sign of the exponent" begin
        K = MortalityTables
        setprecision(BigFloat, 256) do
            logistic(a, k, x) = a * exp(x) / (1 + k * a * exp(x))
            gg(a, b, γ, x) = a * exp(b * x) / (1 + a * γ / b * expm1(b * x))
            # negative growth: exp(-x) overflows where the hazard is still representable
            @test hazard(K.Kannisto(a = 1e308, b = -10.0), 71) ≈ Float64(logistic(big(1e308), 1, big(-710))) rtol = 1e-14
            @test hazard(K.Kannisto(a = 0.5, b = -10.0), 71) ≈ Float64(logistic(big(0.5), 1, big(-710))) rtol = 1e-12
            @test hazard(K.Kannisto(a = 0.5, b = -10.0), 71) > 0   # subnormal, not flushed to zero
            @test hazard(K.Beard(a = 1e308, b = -10.0, k = 0.5), 71) ≈ Float64(logistic(big(1e308), big(0.5), big(-710))) rtol = 1e-14
            @test hazard(K.MakehamBeard(a = 1e308, b = -10.0, k = 0.5), 71) ≈
                  Float64(logistic(big(1e308), big(0.5), big(-710)) + big(0.01)) rtol = 1e-14
            @test hazard(K.KannistoMakeham(a = 1e308, b = -10.0), 71) ≈
                  Float64(logistic(big(1e308), 1, big(-710)) + big(0.001)) rtol = 1e-14
            # positive growth with a large coefficient: k·a would overflow
            @test hazard(K.Beard(a = 1e308, b = 0.1, k = 10.0), 50) ≈ 0.1 rtol = 1e-14
            # GammaGompertz with zero frailty is Gompertz: a decaying hazard, not 0·Inf
            @test hazard(K.GammaGompertz(a = 0.002, b = -10.0, γ = 0.0), 80) == 0.0
            @test hazard(K.GammaGompertz(a = 0.002, b = -10.0, γ = 0.0), Inf) == 0.0
            @test hazard(K.GammaGompertz(a = 0.002, b = -10.0, γ = 0.0), 50) ≈ Float64(big(0.002) * exp(big(-500))) rtol = 1e-12
            @test hazard(K.GammaGompertz(a = 0.5, b = -10.0, γ = 1.0), 71) ≈ Float64(gg(big(0.5), big(-10), big(1), big(71))) rtol = 1e-12
            # Float32 laws keep their type and reach their (subnormal) values
            m32 = K.GammaGompertz(a = 0.002f0, b = -10.0f0, γ = 0.0f0)
            @test hazard(m32, 9.0f0) isa Float32
            @test hazard(m32, 9.0f0) ≈ Float32(0.002 * exp(-90.0)) rtol = 1e-2   # subnormal Float32
            @test hazard(m32, 10.0f0) == 0.0f0
            @test hazard(K.Kannisto(a = 0.5f0, b = -10.0f0), 9.0f0) ≈ Float32(Float64(logistic(big(0.5), 1, big(-90)))) rtol = 1e-2
            @test hazard(K.Kannisto(a = 0.5f0, b = 10.0f0), 100.0f0) == 1.0f0
            # both sides of each branch point: x = 0 for the logistic ratio, and |x| = 1 for
            # GammaGompertz (b = ±1/8 at age 8 is exactly x = ±1)
            for b in (0.125, -0.125), age in (prevfloat(8.0), 8.0, nextfloat(8.0), 0.0, 1e-300)
                @test hazard(K.Beard(a = 0.002, b = b, k = 0.5), age) ≈ Float64(logistic(big(0.002), big(0.5), big(b) * big(age))) rtol = 1e-14
                @test hazard(K.Kannisto(a = 0.5, b = b), age) ≈ Float64(logistic(big(0.5), 1, big(b) * big(age))) rtol = 1e-14
                age > 0 && @test hazard(K.GammaGompertz(a = 0.002, b = b, γ = 1.0), age) ≈
                                 Float64(gg(big(0.002), big(b), big(1), big(age))) rtol = 1e-14
            end
        end
        # Beard with k = 1 is Kannisto: the same hazard, and a quadrature cumulative hazard
        # that agrees with Kannisto's closed form in both tails (it underflowed to a wrong
        # 70.84 for a = 1e308, b = -10 when the hazard was flushed to zero)
        for (a, b) in ((0.5, 10.0), (0.5, -10.0), (1e308, -10.0), (1e308, 0.1)), age in (71.0, 80.0)
            @test hazard(K.Beard(a = a, b = b, k = 1.0), age) == hazard(K.Kannisto(a = a, b = b), age)
            @test cumhazard(K.Beard(a = a, b = b, k = 1.0), age) ≈ cumhazard(K.Kannisto(a = a, b = b), age) rtol = 1e-10
            # the hazard is the derivative of Kannisto's closed-form cumulative hazard
            H(x) = cumhazard(K.Kannisto(a = a, b = b), x)
            h = 1e-5
            @test (H(age + h) - H(age - h)) / 2h ≈ hazard(K.Kannisto(a = a, b = b), age) rtol = 1e-6 atol = 1e-12
        end
        # zero coefficients in Makeham/Gompertz at an infinite age
        @test hazard(Makeham(a = 0.001, b = 0.0, c = 0.002), Inf) == 0.003
        @test hazard(Makeham(a = 0.0, b = 0.13, c = 0.002), Inf) == 0.002
        @test hazard(Gompertz(a = 0.001, b = 0.0), Inf) == 0.001
        @test hazard(Gompertz(a = 0.001, b = -0.1), Inf) == 0.0
        @test hazard(Makeham(), Inf) == Inf
    end

    @testset "parameter derivatives across each algebra switch" begin
        K = MortalityTables
        # a zero amplitude is a smooth boundary at a finite age: the derivatives there are the
        # analytic exp(x) for the hazard and expm1(x)/b for Kannisto's cumulative hazard
        for L in (K.Kannisto, K.Beard, K.MakehamBeard, K.KannistoMakeham, K.GammaGompertz)
            @test ForwardDiff.derivative(a -> hazard(L(a = a, b = 0.13), 50.0), 0.0) ≈ exp(6.5) rtol = 1e-14
            @test ForwardDiff.derivative(a -> hazard(L(a = a, b = -0.13), 50.0), 0.0) ≈ exp(-6.5) rtol = 1e-14
        end
        @test ForwardDiff.derivative(a -> cumhazard(K.Kannisto(a = a, b = 0.13), 50.0), 0.0) ≈ expm1(6.5) / 0.13 rtol = 1e-14
        @test ForwardDiff.derivative(a -> cumhazard(K.KannistoMakeham(a = a, b = 0.13), 50.0), 0.0) ≈ expm1(6.5) / 0.13 rtol = 1e-8
        # the values there are unchanged
        @test hazard(K.Kannisto(a = 0.0, b = 0.13), 50.0) == 0.0
        @test cumhazard(K.Kannisto(a = 0.0, b = 0.13), 50.0) == 0.0

        # Gradients in every parameter against the laws as written, differentiated in 4096-bit
        # arithmetic (the plain formulas cancel at the extreme amplitudes and tiny ages below,
        # so 256 bits isn't enough for a reference). The points sit in and beside each switch
        # of the stable forms: the sign of x = b·age, |x| = 1 (GammaGompertz, and Kannisto's
        # cumulative hazard), a zero amplitude, and amplitudes at which k·a, a·γ/b or a·exp(x)
        # overflows.
        logistic(a, b, k, age) = a * exp(b * age) / (1 + k * a * exp(b * age))
        gg(a, b, γ, age) = a * exp(b * age) / (1 + a * γ / b * expm1(b * age))
        kannisto_H(a, b, age) = log((1 + a * exp(b * age)) / (1 + a)) / b
        function check_gradient(f, ref, p, age)
            ad = ForwardDiff.gradient(q -> f(q, age), p)
            exact = Float64.(ForwardDiff.gradient(q -> ref(q..., big(age)), big.(p)))
            for i in eachindex(p)
                @test ad[i] ≈ exact[i] rtol = 1e-12 atol = 1e-300
            end
        end
        cases = (
            (K.Kannisto, (a, b, age) -> logistic(a, b, 1, age),
                ([0.5, 0.13], [0.0, 0.13], [0.5, -0.13], [0.0, -0.13], [1e308, 0.1])),
            (K.Beard, logistic,   # k·a is 1e308 and then overflows
                ([0.002, 0.13, 0.5], [0.0, 0.13, 0.5], [1e307, 0.1, 10.0], [1e307, 0.1, 20.0])),
            (K.MakehamBeard, (a, b, c, k, age) -> logistic(a, b, k, age) + c,
                ([0.002, 0.13, 0.01, 0.5], [0.0, 0.13, 0.01, 0.5])),
            (K.KannistoMakeham, (a, b, c, age) -> logistic(a, b, 1, age) + c,
                ([0.5, 0.13, 0.001], [0.0, 0.13, 0.001])),
            (K.GammaGompertz, gg,   # a·γ/b is 5e307 and then overflows
                ([0.002, 0.13, 1.0], [0.0, 0.13, 1.0], [0.002, -0.13, 1.0], [1e307, 0.2, 1.0], [1e307, 0.05, 1.0])),
        )
        setprecision(BigFloat, 4096) do
            for (L, ref, params) in cases, p in params, age in (0.0, 5.0, 50.0, 100.0)
                check_gradient((q, age) -> hazard(L(q...), age), ref, p, age)
            end
            # both sides of x = 0 and of |x| = 1 (b = ±1/8 at age 8 is exactly x = ±1)
            for b in (0.125, -0.125), age in (0.0, 1e-300, prevfloat(8.0), 8.0, nextfloat(8.0)), a in (0.0, 0.5)
                check_gradient((q, age) -> hazard(K.Kannisto(q...), age), (a, b, age) -> logistic(a, b, 1, age), [a, b], age)
                check_gradient((q, age) -> hazard(K.Beard(q...), age), logistic, [a, b, 0.5], age)
                check_gradient((q, age) -> hazard(K.GammaGompertz(q...), age), gg, [a, b, 1.0], age)
                age > 0 && check_gradient((q, age) -> cumhazard(K.Kannisto(q...), age), kannisto_H, [a, b], age)
            end
            # Kannisto's cumulative hazard: ordinary and zero amplitudes, negative growth, a large
            # amplitude with negative growth, and each side of a·exp(x) overflowing (x = 709, 710)
            for (p, age) in (
                    ([0.5, 0.13], 5.0), ([0.5, 0.13], 50.0), ([0.5, 0.13], 100.0), ([0.0, 0.13], 50.0),
                    ([0.5, -0.13], 50.0), ([0.0, -0.13], 50.0), ([1e308, -10.0], 71.0),
                    ([0.5, 10.0], 70.9), ([0.5, 10.0], 71.0),
                )
                check_gradient((q, age) -> cumhazard(K.Kannisto(q...), age), kannisto_H, p, age)
            end
        end
    end

    @testset "Makeham" begin

        g = Gompertz(a=0.0002, b=.13)

        @test survival(g, 45) ≈ 0.5870365016720939
        @test survival(g, 45, 46) ≈ 0.9285202788707242
        @test decrement(g, 45, 46) ≈ 1 - 0.9285202788707242
        @test hazard(g, 45) ≈ 0.06944687609574696

        @testset "AMLCR" begin
            m = Makeham(a=2.7e-6, b=log(1.124), c=0.00022)

            @test MortalityTables.μ(m, 20) == 0.00022 + 2.7e-6 * 1.124^20
            @test_throws MethodError m[20]  # indexing a model is not supported; call it or use hazard
            @test m(20) == MortalityTables.μ(m, 20)
            @test MortalityTables.μ === hazard
            
            # vs manually calculated (via QuadGK) integrals
            @test decrement(m, 20, 25) ≈ 0.0012891622754368504
            @test survival(m, 20, 25) ≈ 1 - 0.0012891622754368504
            @test decrement(m, 25) ≈ 0.005888764668801838
            @test survival(m, 25) ≈ 1 - 0.005888764668801838

            # these values come from the 'Standard Select and Ultimate Survival Model'
            # from Actuarial Mathematics for Life Contingent Risks, 2nd end

            ℓ = 100_000
            ℓs = [survival(m, 20, age) for age in 21:100] .* ℓ

            ℓ_age(x) = ℓs[x - 20]
            @test isapprox(ℓ_age(21),  99975.04, atol=0.01)
            @test isapprox(ℓ_age(31),  99695.83, atol=0.01)
            @test isapprox(ℓ_age(82),  70507.19, atol=0.01)
            @test isapprox(ℓ_age(100),  6248.17, atol=0.01)
        end

    end

    @testset "table stand-in semantics" begin
        m = Makeham()
        # a DeathDistribution is accepted and ignored by continuous models
        @test survival(m, 65, Uniform()) == survival(m, 65)
        @test survival(m, 60, 65, Uniform()) == survival(m, 60, 65)
        @test decrement(m, 60, 65) ≈ 1 - survival(m, 60, 65)
        @test decrement(m, 65, Uniform()) == decrement(m, 65)
        @test decrement(m, 60, 65, Uniform()) == decrement(m, 60, 65)
        # a small decrement is not rounded to zero by 1 - survival
        tiny = Makeham(a = 1e-12, b = 0.1, c = 0.0)
        @test decrement(tiny, 1e-6) ≈ cumhazard(tiny, 1e-6) rtol = 1e-12
        @test decrement(tiny, 1e-6) > 0
        @test omega(m) == Inf
        # a law whose formula ends has a finite omega, with survival still positive there
        w = MortalityTables.Wittstein()
        @test omega(w) == 100
        @test isfinite(hazard(w, 100)) && survival(w, 100) > 0
        @test_throws DomainError hazard(w, 101)
        @test omega(MortalityTables.VanderMaen()) == 200
        @test omega(MortalityTables.VanderMaen2(n = 150)) == 150
        @test omega(MortalityTables.Wittstein(m = 90.0)) === 90.0
    end

    @testset "Gompertz and Makeham equality" begin

        # Gompertz is Makeham's where c = 0
        m = Makeham(a=2.7e-6, b=1.124, c=0.0)
        g = Gompertz(a=2.7e-6, b=1.124)

        for age ∈ 20:100
            @test survival(m, age) == survival(g, age)
            @test survival(m, age, 1) == survival(g, age, 1)
        end
    end

    @testset "MortalityLaws R package" begin
        model_tests = [(rmodel = "gompertz", juliamodel = MortalityTables.Gompertz()),
                        # (rmodel="gompertz0",juliamodel=MortalityTables.nothing),
                        (rmodel = "invgompertz", juliamodel = MortalityTables.InverseGompertz()),
                        (rmodel = "makeham", juliamodel = MortalityTables.Makeham()),
                        # (rmodel="makeham0",juliamodel=MortalityTables.nothing),
                        (rmodel = "opperman", juliamodel = MortalityTables.Opperman()),
                        (rmodel = "thiele", juliamodel = MortalityTables.Thiele()),
                        (rmodel = "wittstein", juliamodel = MortalityTables.Wittstein()),
                        (rmodel = "perks", juliamodel = MortalityTables.Perks()),
                        (rmodel = "weibull", juliamodel = MortalityTables.Weibull()),
                        (rmodel = "invweibull", juliamodel = MortalityTables.InverseWeibull()),
                        (rmodel = "vandermaen", juliamodel = MortalityTables.VanderMaen()),
                        (rmodel = "vandermaen2", juliamodel = MortalityTables.VanderMaen2()),
                        (rmodel = "strehler_mildvan", juliamodel = MortalityTables.StrehlerMildvan()),
                        (rmodel = "quadratic", juliamodel = MortalityTables.Quadratic()),
                        (rmodel = "beard", juliamodel = MortalityTables.Beard()),
                        (rmodel = "beard_makeham", juliamodel = MortalityTables.MakehamBeard()),
                        (rmodel = "ggompertz", juliamodel = MortalityTables.GammaGompertz()),
                        (rmodel = "siler", juliamodel = MortalityTables.Siler()),
                        (rmodel = "HP", juliamodel = MortalityTables.HeligmanPollard()),
                        (rmodel = "HP2", juliamodel = MortalityTables.HeligmanPollard2()),
                        (rmodel = "HP3", juliamodel = MortalityTables.HeligmanPollard3()),
                        (rmodel = "HP4", juliamodel = MortalityTables.HeligmanPollard4()),
                        (rmodel = "rogersplanck", juliamodel = MortalityTables.RogersPlanck()),
                        (rmodel = "martinelle", juliamodel = MortalityTables.Martinelle()),
                        (rmodel = "kostaki", juliamodel = MortalityTables.Kostaki()),
                        (rmodel = "kannisto", juliamodel = MortalityTables.Kannisto()),
                        (rmodel = "kannisto_makeham", juliamodel = MortalityTables.KannistoMakeham())
                        # the next two requre adding an autodiff dependency:
                        # (rmodel="carriere1",juliamodel=MortalityTables.Carriere()),
                        # (rmodel="carriere2",juliamodel=MortalityTables.Carriere2()),
                    ]


        # load test targets from data
        dir = joinpath(pwd(), "data", "parametric")

        rmodels = Dict()
        for (root, dirs, files) in walkdir(dir)
            for file in files
                if file[end - 3:end] == "json"
                    json = JSON.parse(MortalityTables.open_and_read(joinpath(root, file)))
                    rmodels[json["modelname"][1]] = json
                end
            end
        end

        # compare values
        ages = 20:100
        @testset "Model: $(model.rmodel)" for model in model_tests
            rmodel = rmodels[model.rmodel]
            # test hazard
            if "hx" in keys(rmodel)
                for (i, age) in enumerate(ages)
                    @test hazard(model.juliamodel, age) ≈ rmodel["hx"][i] 
                end
            end

            # test cumulative hazard

            if "Hx" in keys(rmodel)
                for (i, age) in enumerate(ages)
                    @test cumhazard(model.juliamodel, age) ≈ rmodel["Hx"][i] 
                end
            end
            
            # test Survival
            if "Sx" in keys(rmodel)
                for (i, age) in enumerate(ages)
                    @test survival(model.juliamodel, age) ≈ rmodel["Sx"][i] 
                end
            end

            # test other characteristics

            @test survival(model.juliamodel, 20, 20) == 1.0
            if model.rmodel == "perks"
                # the default params give a hazard that is numerically zero at these ages
                @test_broken survival(model.juliamodel, 50, 51) < 1.0
            else
                @test survival(model.juliamodel, 50, 51) < 1.0
            end
            
            @test model.juliamodel(20) >= 0

        end

    end

end