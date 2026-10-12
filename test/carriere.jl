using ForwardDiff

@testset "Alternative laws and Carrière mixtures" begin
    K = MortalityTables

    @testset "location/scale parameterization" begin
        for (μ, σ, c) in ((49.0, 7.7, 0.001), (20.0, 3.0, 0.0), (-5.0, 8.0, 0.1))
            m = K.Makeham2(μ = μ, σ = σ, c = c)
            ref = Makeham(a = exp(-μ / σ) / σ, b = 1 / σ, c = c)
            for age in (0.0, 1e-14, 0.5, 20.0, 80.0, 150.0)
                @test hazard(m, age) ≈ hazard(ref, age) rtol = 2e-14
                @test cumhazard(m, age) ≈ cumhazard(ref, age) rtol = 2e-14
            end
            @test cumhazard(m, 0) == 0
            @test survival(m, 0) == 1
            @test cumhazard(m, 50, 50) == 0
            @test cumhazard(m, 50, 40) ≈ -cumhazard(m, 40, 50)
            @test ForwardDiff.derivative(x -> cumhazard(m, x), 0.0) ≈ hazard(m, 0)
        end
        @test K.Gompertz2(μ = 40, σ = 8) === K.Makeham2(μ = 40, σ = 8, c = 0)
        @test K.Gompertz2(mu = 40.0f0, sigma = 8.0f0) isa K.Makeham2{Float32}
        @test K.Makeham2(mu = 40, sigma = 8) === K.Makeham2(μ = 40, σ = 8)
        for L in (K.Makeham2, K.Gompertz2)
            @test_throws ArgumentError L(μ = 40, mu = 40)
            @test_throws ArgumentError L(σ = 8, sigma = 8)
            for bad in (0.0, -1.0, Inf, NaN)
                @test_throws ArgumentError L(σ = bad)
            end
            for bad in (Inf, -Inf, NaN)
                @test_throws ArgumentError L(μ = bad)
            end
        end
        for bad in (-1.0, Inf, NaN)
            @test_throws ArgumentError K.Makeham2(c = bad)
        end
        @test_throws ArgumentError K.Makeham2(49, -1, 0)

        # The equivalent Makeham coefficient exp(-μ/σ)/σ would underflow to zero.
        m = K.Gompertz2(μ = 1000, σ = 1)
        @test hazard(m, 1000) ≈ 1
        @test cumhazard(m, 1000) ≈ 1
        @test survival(m, 1000) ≈ exp(-1)
        @test hazard(m, 1100) ≈ exp(100)
        @test cumhazard(m, 1100) ≈ exp(100)
        @test survival(m, 1100, 1101) == 0
        @test survival(m, 2000, 2000) == 1
        for m in (K.Makeham2(), K.Gompertz2())
            @test survival(m, Inf) == 0
            @test hazard(m, Inf) == Inf
        end
    end

    @testset "independent MortalityLaws component fixtures" begin
        # These stored R fixtures predate this implementation. Combine the independently
        # computed component S and h values, not carriere*.json: those contain finite
        # differences, and carriere1 also uses sigma1 for its inverse-Weibull component.
        # The explicit parameters below match each fixture's par/param_names metadata.
        fixture(name) = JSON.parsefile(joinpath(@__DIR__, "data", "parametric", name * ".json"))
        wei, iw, ig, gom = map(fixture, ("weibull", "invweibull", "invgompertz", "gompertz0"))
        cases = ((K.Carriere1(μ₁ = 1, σ₁ = 2, μ₂ = 5, σ₂ = 10, μ₃ = 49, σ₃ = 7.7), (wei, iw, gom)),
                 (K.Carriere2(σ₂ = 7.7, σ₃ = 7.7), (wei, ig, gom)))
        for (m, fixtures) in cases, (i, age) in enumerate(20:100)
            weights = (m.p₁, m.p₂, 1 - m.p₁ - m.p₂)
            S = sum(w * f["Sx"][i] for (w, f) in zip(weights, fixtures))
            h = sum(w * f["Sx"][i] * f["hx"][i] for (w, f) in zip(weights, fixtures)) / S
            @test survival(m, age) ≈ S rtol = 2e-13
            @test cumhazard(m, age) ≈ -log(S) rtol = 2e-13
            @test hazard(m, age) ≈ h rtol = 2e-13
        end
    end

    @testset "calculus and interval contract" begin
        # Non-default scales catch a component using another component's parameter.
        for m in (K.Carriere1(), K.Carriere2(),
                  K.Carriere1(p₁ = 0.2, p₂ = 0.3, μ₁ = 2, σ₁ = 1, μ₂ = 8, σ₂ = 3),
                  K.Carriere2(p₁ = 0.2, p₂ = 0.3, μ₁ = 2, σ₁ = 1, μ₂ = 20, σ₂ = 3))
            @test cumhazard(m, 0) == 0
            @test survival(m, 0) == 1
            for age in (0.1, 2.0, 20.0, 60.0, 100.0)
                δ = 1e-5 * age
                fd = (cumhazard(m, age + δ) - cumhazard(m, age - δ)) / (2δ)
                @test fd ≈ hazard(m, age) rtol = 1e-7
                @test ForwardDiff.derivative(x -> cumhazard(m, x), age) ≈ hazard(m, age) rtol = 1e-11
                @test cumhazard(m, age) ≈ quadgk(x -> hazard(m, x), 0, age; rtol = 1e-9)[1] rtol = 1e-8
            end
            ages = [0.0, 0.1, 3.5, 20.0, 81.0, 110.0]
            @test hazard.(m, ages) == [hazard(m, x) for x in ages]
            @test issorted(cumhazard.(m, ages))
            @test issorted(survival.(m, ages); rev = true)
            @test all(x -> 0 <= x <= 1, survival.(m, ages))
            @test survival(m, 20, 40) ≈ survival(m, 40) / survival(m, 20)
            @test survival(m, 20, 80) ≈ survival(m, 20, 40) * survival(m, 40, 80)
            @test survival(m, 80, 20) ≈ inv(survival(m, 20, 80))
            @test survival(m, 80, 80) == 1
            @test decrement(m, 20, 80) ≈ 1 - survival(m, 20, 80)
            @test survival(m, 20, 80, ConstantForce()) == survival(m, 20, 80)
            @test m(20) == hazard(m, 20)
        end
        @test survival(K.Carriere1(σ₂ = 6), 20) != survival(K.Carriere1(σ₂ = 15), 20)
        @test survival(K.Carriere2(σ₂ = 6), 20) != survival(K.Carriere2(σ₂ = 15), 20)
        for L in (K.Carriere1, K.Carriere2)
            m = L()
            args = collect(map(f -> getfield(m, f), fieldnames(typeof(m))))
            f(v) = cumhazard(L(v...), 40.0)
            gradient = ForwardDiff.gradient(f, args)
            @test all(isfinite, gradient)
            for i in eachindex(args)
                δ = args[i] * 1e-5
                lo, hi = copy(args), copy(args)
                lo[i] -= δ
                hi[i] += δ
                @test gradient[i] ≈ (f(hi) - f(lo)) / (2δ) rtol = 1e-5 atol = 1e-12
            end
        end
    end

    @testset "simplex boundaries and zero-age limits" begin
        for L in (K.Carriere1, K.Carriere2)
            m = L()
            middle = L === K.Carriere1 ? K.InverseWeibull(m = m.μ₂, σ = m.σ₂) : K.InverseGompertz(m = m.μ₂, σ = m.σ₂)
            cases = ((L(p₁ = 1, p₂ = 0), K.Weibull(m = m.μ₁, σ = m.σ₁)),
                     (L(p₁ = 0, p₂ = 1), middle),
                     (L(p₁ = 0, p₂ = 0), K.Gompertz2(μ = m.μ₃, σ = m.σ₃)))
            for (mixture, component) in cases, age in (0.0, 1e-12, 20.0, 100.0, 1e4, Inf)
                @test cumhazard(mixture, age) ≈ cumhazard(component, age)
                @test hazard(mixture, age) ≈ hazard(component, age)
            end
            for (p1, p2) in ((0.0, 0.4), (0.2, 0.0), (0.2, 0.8), (0.0, 0.0))
                model = L(p₁ = p1, p₂ = p2)
                @test survival(model, 0) == 1
                @test !isnan(hazard(model, 0))
                @test !isnan(hazard(model, 1e4))
            end
            @test hazard(L(), 0) == Inf
            @test hazard(L(p₁ = 1e-20, p₂ = 0), 0) == Inf
            @test isfinite(hazard(L(p₁ = 0), 0))
            # Weibull's finite/non-singular limits, with the other weights zero.
            @test hazard(L(p₁ = 1, p₂ = 0, μ₁ = 2, σ₁ = 2), 0) ≈ 0.5
            @test hazard(L(p₁ = 1, p₂ = 0, μ₁ = 3, σ₁ = 2), 0) == 0
            @test cumhazard(L(), Inf) == Inf
            @test survival(L(), Inf) == 0
            @test hazard(L(), Inf) == 0
            @test hazard(L(p₁ = 0, p₂ = 0), Inf) == Inf
            @test survival(L(p₁ = 0, p₂ = 0), 10000, 10001) == 0
            @test survival(L(p₁ = 0, p₂ = 0), 10000, 10000) == 1
        end
        @test hazard(K.Carriere2(p₁ = 0, p₂ = 0.5), Inf) ≈ 1 / 7
        @test hazard(K.Carriere2(p₁ = 0.3, μ₁ = 2, σ₁ = 1), Inf) ≈ 1 / 7
    end

    @testset "validation and aliases" begin
        for L in (K.Carriere1, K.Carriere2)
            for (p1, p2) in ((-0.01, 0.2), (0.2, -0.01), (0.8, 0.8), (Inf, 0.0), (0.0, NaN))
                @test_throws ArgumentError L(p₁ = p1, p₂ = p2)
                @test_throws ArgumentError L(p1, p2, 1, 2, 3, 4, 5, 6)
            end
            for field in (:μ₁, :μ₂, :μ₃, :σ₁, :σ₂, :σ₃), bad in (0.0, -1.0, Inf, NaN)
                @test_throws ArgumentError L(; field => bad)
            end
            for (unicode, ascii, value) in ((:p₁, :p1, 0.1), (:p₂, :p2, 0.2),
                (:μ₁, :mu1, 2), (:μ₂, :mu2, 4), (:μ₃, :mu3, 90),
                (:σ₁, :sigma1, 3), (:σ₂, :sigma2, 7), (:σ₃, :sigma3, 10))
                @test L(; unicode => value) === L(; ascii => value)
                @test_throws ArgumentError L(; unicode => value, ascii => value)
            end
        end
        for m in (K.Makeham2(), K.Gompertz2(), K.Carriere1(), K.Carriere2()), age in (-1.0, NaN)
            @test_throws DomainError hazard(m, age)
            @test_throws DomainError cumhazard(m, age)
            @test_throws DomainError survival(m, age)
            @test_throws DomainError decrement(m, 0, age)
        end
    end

    @testset "small decrements and tails against high precision" begin
        # Independent, direct survival formulas in 1024-bit arithmetic. Production
        # evaluates cumulative hazards and scaled masses instead of these expressions.
        setprecision(BigFloat, 1024) do
            function reference(m, x)
                p1, p2, m1, m2, m3, s1, s2, s3 =
                    map(f -> BigFloat(getfield(m, f)), fieldnames(typeof(m)))
                x = BigFloat(x)
                w = (p1, p2, 1 - p1 - p2)
                S1 = exp(-(x / m1)^(m1 / s1))
                S2 = m isa K.Carriere1 ? -expm1(-(x / m2)^(-m2 / s2)) :
                     -expm1(-exp((m2 - x) / s2)) / -expm1(-exp(m2 / s2))
                S3 = exp(-exp(-m3 / s3) * expm1(x / s3))
                h1 = (x / m1)^(m1 / s1 - 1) / s1
                h2 = m isa K.Carriere1 ? (x / m2)^(-m2 / s2 - 1) / (s2 * expm1((x / m2)^(-m2 / s2))) :
                     exp((m2 - x) / s2) / (s2 * expm1(exp((m2 - x) / s2)))
                h3 = exp((x - m3) / s3) / s3
                S = w[1] * S1 + w[2] * S2 + w[3] * S3
                h = (w[1] * S1 * h1 + w[2] * S2 * h2 + w[3] * S3 * h3) / S
                return -log(S), h
            end
            for m in (K.Carriere1(), K.Carriere2(), K.Carriere1(p₁ = 0), K.Carriere2(p₁ = 0)),
                age in (1e-12, 0.1, 100.0, 1000.0, 10000.0)
                H, h = reference(m, age)
                @test cumhazard(m, age) ≈ H rtol = 2e-12 atol = 1e-300
                @test hazard(m, age) ≈ h rtol = 2e-12 atol = 1e-300
            end
            for m in (K.Carriere1(), K.Carriere2(), K.Carriere2(p₁ = 0))
                a, b = 10000.0, 10001.0
                Ha, _ = reference(m, a)
                Hb, _ = reference(m, b)
                @test survival(m, a, b) ≈ exp(Ha - Hb) rtol = 1e-12
                @test decrement(m, a, b) ≈ -expm1(Ha - Hb) rtol = 2e-10
            end
        end
        # The inverse laws must retain their own tiny hazards and distant tails too.
        for m in (K.InverseGompertz(m = 1, σ = 2), K.InverseWeibull(m = 1, σ = 2))
            @test cumhazard(m, 0) == 0
            @test survival(m, 0) == 1
            @test isfinite(hazard(m, 1e4))
            @test isfinite(cumhazard(m, 1e4))
            @test survival(m, Inf) == 0
        end
        m = K.InverseGompertz(m = 1, σ = 2)
        @test cumhazard(m, 1e-20) ≈ hazard(m, 0) * 1e-20 rtol = 1e-14
        @test hazard(K.InverseWeibull(), 0) == 0
        @test hazard(K.InverseWeibull(), Inf) == 0
        @test cumhazard(K.InverseWeibull(m = 1, σ = 2), 1e300) ≈ 0.5 * log(1e300)
        m = K.Carriere2(p₁ = 0, p₂ = 1)
        @test survival(m, 1e4) == 0
        @test survival(m, 1e4, 1e4 + 1) ≈ exp(-1 / 7)
        @test decrement(K.Carriere1(p₁ = 0, p₂ = 0), 1e-20) > 0
        # All active cumulative hazards overflow, but the mixture hazard and a
        # sufficiently short interval hazard can still be represented. The dominant
        # component switches between these ages; compare H, not instantaneous h.
        m = K.Carriere2(p₁ = 0.4, p₂ = 0.6, μ₁ = 1, σ₁ = 0.5, μ₂ = 1, σ₂ = 1e-156)
        for (age, h) in ((1e155, 2e155), (1e157, 1e156))
            @test cumhazard(m, age) == Inf
            @test hazard(m, age) ≈ h rtol = 1e-12
            to = nextfloat(age)
            H = age == 1e155 ? BigFloat(to)^2 - BigFloat(age)^2 : BigFloat(to - age) / BigFloat(m.σ₂)
            @test cumhazard(m, age, to) ≈ H rtol = 1e-12
            @test survival(m, age, to) == 0
        end
        # An inverse-Gompertz law approaches an exponential as its location tends
        # far below zero; normalization must not become 0/0 in that limit.
        m = K.InverseGompertz(m = -1000, σ = 1)
        @test cumhazard(m, 1e-20) ≈ 1e-20
        @test hazard(m, 0) ≈ 1
    end
end
