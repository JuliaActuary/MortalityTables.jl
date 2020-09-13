@testset "Parameterized Models" begin

    @testset "Makeham" begin

        g = Gompertz(a=0.0002,b=.13)

        @test survivorship(g,45) ≈ 0.5870365016720939
        @test survivorship(g,45,46) ≈ 0.9285202788707242
        @test decrement(g,45,46) ≈ 1 - 0.9285202788707242
        @test hazard(g,45) ≈ 0.06944687609574696

        @testset "AMLCR" begin
            m = Makeham(a=2.7e-6,b=log(1.124), c= 0.00022)

            @test MortalityTables.μ(m, 20) == 0.00022 + 2.7e-6 * 1.124^20
            @test m[20] == MortalityTables.μ(m, 20)
            @test m(20) == MortalityTables.μ(m, 20)
            
            # vs manually calculated (via QuadGK) integrals
            @test decrement(m, 20, 25) ≈ 0.0012891622754368504
            @test survivorship(m, 20, 25) ≈ 1 - 0.0012891622754368504
            @test decrement(m, 25) ≈ 0.005888764668801838
            @test survivorship(m, 25) ≈ 1 - 0.005888764668801838

            # these values come from the 'Standard Select and Ultimate Survival Model'
            # from Actuarial Mathematics for Life Contingent Risks, 2nd end

            ℓ = 100_000
            ℓs = [survivorship(m, 20, age) for age in 21:100] .* ℓ

            ℓ_age(x) = ℓs[x - 20]
            @test isapprox(ℓ_age(21),  99975.04, atol = 0.01)
            @test isapprox(ℓ_age(31),  99695.83, atol = 0.01)
            @test isapprox(ℓ_age(82),  70507.19, atol = 0.01)
            @test isapprox(ℓ_age(100),  6248.17, atol = 0.01)
        end

    end

    @testset "Gompertz and Makeham equality" begin

        # Gompertz is Makeham's where c = 0
        m = Makeham( a=2.7e-6,b= 1.124,c=0.0)
        g = Gompertz(a=2.7e-6,b= 1.124)

        for age ∈ 20:100
            @test survivorship(m, age) == survivorship(g, age)
            @test survivorship(m, age, 1) == survivorship(g, age, 1)
        end
    end

    @testset "MortalityLaws R package" begin
        model_tests = [(rmodel="gompertz",juliamodel=MortalityTables.Gompertz(),test=true),
                        (rmodel="gompertz0",juliamodel=MortalityTables.nothing,test=false),
                        (rmodel="invgompertz",juliamodel=MortalityTables.InverseGompertz(),test=true),
                        (rmodel="makeham",juliamodel=MortalityTables.Makeham(),test=true),
                        (rmodel="makeham0",juliamodel=MortalityTables.nothing,test=false),
                        (rmodel="opperman",juliamodel=MortalityTables.Opperman(),test=true),
                        (rmodel="thiele",juliamodel=MortalityTables.Thiele(),test=true),
                        (rmodel="wittstein",juliamodel=MortalityTables.Wittstein(),test=true),
                        (rmodel="perks",juliamodel=MortalityTables.Perks(),test=true),
                        (rmodel="weibull",juliamodel=MortalityTables.Weibull(),test=true),
                        (rmodel="invweibull",juliamodel=MortalityTables.InverseWeibull(),test=true),
                        (rmodel="vandermaen",juliamodel=MortalityTables.VanderMaen(),test=true),
                        (rmodel="vandermaen2",juliamodel=MortalityTables.VanderMaen2(),test=true),
                        (rmodel="strehler_mildvan",juliamodel=MortalityTables.StrehlerMildvan(),test=true),
                        (rmodel="quadratic",juliamodel=MortalityTables.Quadratic(),test=true),
                        (rmodel="beard",juliamodel=MortalityTables.Beard(),test=true),
                        (rmodel="beard_makeham",juliamodel=MortalityTables.MakehamBeard(),test=true),
                        (rmodel="ggompertz",juliamodel=MortalityTables.GammaGompertz(),test=false),
                        (rmodel="siler",juliamodel=MortalityTables.nothing,test=false),
                        (rmodel="HP",juliamodel=MortalityTables.HeligmanPollard(),test=false),
                        (rmodel="HP2",juliamodel=MortalityTables.HeligmanPollard2(),test=false),
                        (rmodel="HP3",juliamodel=MortalityTables.HeligmanPollard3(),test=false),
                        (rmodel="HP4",juliamodel=MortalityTables.HeligmanPollard4(),test=false),
                        (rmodel="rogersplanck",juliamodel=MortalityTables.RogersPlanck(),test=false),
                        # (rmodel="martinelle",juliamodel=MortalityTables.Martinelle(),test=false),
                        # (rmodel="kostaki",juliamodel=MortalityTables.Kostaki(),test=false),
                        # (rmodel="carriere1",juliamodel=MortalityTables.Carriere(),test=false),
                        # (rmodel="carriere2",juliamodel=MortalityTables.Carriere2(),test=false),
                        # (rmodel="kannisto",juliamodel=MortalityTables.Kannisto(),test=false),
                        # (rmodel="kannisto_makeham",juliamodel=MortalityTables.KannistoMakeham(),test=false)
                    ]


        # load test targets from data
        dir = joinpath(pwd(),"data","parametric")

        rmodels = Dict()
        for (root, dirs, files) in walkdir(dir)
            for file in files
                if file[end-3:end] == "json"
                    json = JSON.parse(MortalityTables.open_and_read(joinpath(root, file)))
                    rmodels[json["modelname"][1]] = json
                end
            end
        end

        # compare values
        ages = 20:100
        @testset "Model: $(model.rmodel)" for model in model_tests
            if model.test || ~isnothing(model.juliamodel)
                rmodel = rmodels[model.rmodel]
                # test hazard
                if "hx" in keys(rmodel)
                    for (i,age) in enumerate(ages)
                        @test hazard(model.juliamodel,age) ≈ rmodel["hx"][i] 
                    end
                end

                # test cumulative hazard

                if "Hx" in keys(rmodel)
                    for (i,age) in enumerate(ages)
                        @test cumhazard(model.juliamodel,age) ≈ rmodel["Hx"][i] 
                    end
                end
                
                #test Survival
                if "Sx" in keys(rmodel)
                    for (i,age) in enumerate(ages)
                        @test survivorship(model.juliamodel,age) ≈ rmodel["Sx"][i] 
                    end
                end


            end
        end

    end

end