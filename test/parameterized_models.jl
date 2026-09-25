using InteractiveUtils: subtypes

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
            @test decrement(m, 40, 60) ≈ 1 - survival(m, 40, 60)
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
        @test decrement(m, 60, 65) == 1 - survival(m, 60, 65)
        @test decrement(m, 65, Uniform()) == 1 - survival(m, 65)
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