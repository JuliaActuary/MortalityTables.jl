@testset "MortalityTables from File" begin

    @testset "TableMetaData" begin
        @testset "whitespace management" begin
            # "2017 Loaded CSO Smoker Distinct Nonsmoker Female ALB" comes with a trailing whitespace
            @test "2017 Loaded CSO Smoker Distinct Nonsmoker Female ALB" in
                  keys(tables)
        end
    end

    @testset "SOA tables" begin

        @testset "Ultimate Only SOA table" begin
            cso1980 = tables["1980 CSO Basic Table – Male, ANB"]
            @test q(cso1980.ultimate, 35, 1) ≈ 0.00118
            @test q(cso1980, 35, 1) ≈ 0.00118
            @test q(cso1980.ultimate, 35, 1:2) ≈ [0.00118, 0.00128]
            @test q(cso1980.ultimate, 35, 1:2) ≈ [0.00118, 0.00128]
            @test q(cso1980.ultimate, 35, 1) ≈ 0.00118
            @test q(cso1980.ultimate, 35, 61) ≈ 0.27302
            @test q(cso1980.ultimate, 95, 1) ≈ 0.27302
            @test_throws BoundsError q(cso1980.ultimate, 35, 95)
            @test ismissing(q(cso1980.ultimate, 101, 1))
            @test omega(cso1980.ultimate, 35) == 100
            @test ω(cso1980.ultimate, 35) == 100
            @test ω(cso1980, 35) == 100
        end

        @testset "2001 VBT" begin
            vbt2001 =
                tables["2001 VBT Residual Standard Select and Ultimate - Male Nonsmoker, ANB"]


            @test q(vbt2001.select, 35, 1) ≈ 0.00036
            @test q(vbt2001.select, 35, 61) ≈ 0.24298
            @test q(vbt2001.select, 95, 1) ≈ 0.23815
            @test q(vbt2001.ultimate, 95) ≈ 0.24298
            @test q(vbt2001.ultimate, 120) ≈ 1.0
            @test q(vbt2001.select, 120, 1) ≈ 1.0
            @test_throws BoundsError q(vbt2001.select, 35, 95)
            @test ismissing(q(vbt2001.select, 150, 1))
            @test ismissing(q(vbt2001.ultimate, 150, 1))
            @test omega(vbt2001.select, 20) == 120
            @test ismissing(ω(vbt2001.ultimate, 20))
            @test ω(vbt2001.select, 25) == 120

            # this tests that mortality vectors are generated
            # when there are ultimate rates available
            # but no select rates for the issue age
            @test q(vbt2001.select, 116, 1) ≈ 0.79988

        end

        @testset "2001 CSO" begin
            cso2001 =
                tables["2001 CSO Super Preferred Select and Ultimate - Male Nonsmoker, ANB"]

            @test q(cso2001.select, 35, 1) ≈ 0.00037
            @test q(cso2001.select, 35, 61) ≈ 0.26719
            @test q(cso2001.ultimate, 16, 1) ≈ 0.00041
            @test q(cso2001.ultimate, 95) ≈ 0.26719
            @test q(cso2001.ultimate, 95, 1) ≈ 0.26719
            @test q(cso2001.ultimate, 120) ≈ 1.0
            @test q(cso2001.ultimate, 120, 1) ≈ 1.0
            @test q(cso2001.select, 120, 1) ≈ 1.0
            @test ismissing(q(cso2001.select, 15, 1))
            @test ismissing(q(cso2001.ultimate, 150))
            @test_throws BoundsError q(cso2001.select, 35, 95)
            @test omega(cso2001.select, 16) == 120
            @test ω(cso2001.ultimate, 16) == 120
            @test ismissing(ω(cso2001.ultimate, 10))

            @testset "duration ranges" begin
                #select only
                @test q(cso2001.select, 29, 1:3) == [0.00029, 0.00035, 0.0004]
                @test q(cso2001.select, 35, 1:25) == [
                    0.00037,
                    0.00043,
                    0.00049,
                    0.00057,
                    0.00063,
                    0.0007,
                    0.00077,
                    0.00084,
                    0.00092,
                    0.00101,
                    0.00114,
                    0.00127,
                    0.00143,
                    0.00159,
                    0.00174,
                    0.00188,
                    0.00208,
                    0.00231,
                    0.00251,
                    0.00279,
                    0.00315,
                    0.00357,
                    0.00406,
                    0.00461,
                    0.00508,
                ]

                # crosses into ultimate
                @test q(cso2001.select, 35, 21:30) == [
                    0.00315,
                    0.00357,
                    0.00406,
                    0.00461,
                    0.00508,
                    0.00621,
                    0.0069,
                    0.00773,
                    0.00867,
                    0.00965,
                ]
                @test q(cso2001.select, 35, 1:30) == [
                    0.00037,
                    0.00043,
                    0.00049,
                    0.00057,
                    0.00063,
                    0.0007,
                    0.00077,
                    0.00084,
                    0.00092,
                    0.00101,
                    0.00114,
                    0.00127,
                    0.00143,
                    0.00159,
                    0.00174,
                    0.00188,
                    0.00208,
                    0.00231,
                    0.00251,
                    0.00279,
                    0.00315,
                    0.00357,
                    0.00406,
                    0.00461,
                    0.00508,
                    0.00621,
                    0.0069,
                    0.00773,
                    0.00867,
                    0.00965,
                ]
            end

            @testset "Time Based Calcs" begin
                @test q(cso2001.select, 29, 1, 1) == 1.0 - prod(1 .- [0.00029])
                @test q(cso2001.select, 29, 1, 2) ==
                      1.0 - prod(1 .- [0.00029, 0.00035])
                @test p(cso2001.select, 29, 1, 2) ==
                      prod(1 .- [0.00029, 0.00035])

            end
        end





    end
end
