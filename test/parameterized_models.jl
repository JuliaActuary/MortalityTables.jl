@testset "Parameterized Models" begin

    @testset "Makeham" begin

        g = Gompertz(a=0.0002, b=.13)

        @test survival(g, 45) ≈ 0.5870365016720939
        @test survival(g, 45, 46) ≈ 0.9285202788707242
        @test decrement(g, 45, 46) ≈ 1 - 0.9285202788707242
        @test hazard(g, 45) ≈ 0.06944687609574696

        @testset "AMLCR" begin
            m = Makeham(a=2.7e-6, b=log(1.124), c=0.00022)

            @test MortalityTables.μ(m, 20) == 0.00022 + 2.7e-6 * 1.124^20
            @test m[20] == MortalityTables.μ(m, 20)
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
            if model.rmodel in ["quadratic","perks","vandermaen","vandermaen2"]
                # the default params create a crazy hazard function 
                @test_broken survival(model.juliamodel, 50, 51) < 1.0
            else
                @test survival(model.juliamodel, 50, 51) < 1.0
            end
            
            @test model.juliamodel[20] >= 0

        end

    end

end