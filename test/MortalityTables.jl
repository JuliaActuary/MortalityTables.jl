@testset "MortalityTables from File" begin

    @testset "TableMetaData" begin
        @testset "whitespace management" begin
            # "2017 Loaded CSO Smoker Distinct Nonsmoker Female ALB" comes with a trailing whitespace
            @test ~isnothing(MortalityTables.table_source_map["2017 Loaded CSO Smoker Distinct Nonsmoker Female ALB"])
        end
    end

    @testset "SOA tables" begin

        @testset "Ultimate Only SOA table" begin
            cso1980 = MortalityTables.table("1980 CSO Basic Table – Male, ANB")
            @test cso1980.ultimate[35] ≈ 0.00118
            @test cso1980[35] ≈ 0.00118
            @test cso1980[35:36] ≈ [0.00118, 0.00128]
            @test cso1980[35 + 60] ≈ 0.27302
            @test cso1980[95] ≈ 0.27302
            @test_throws BoundsError cso1980[125]
            @test_throws ArgumentError survival(cso1980,10,15)
            @test_throws ArgumentError decrement(cso1980,10,15)
            @test omega(cso1980.ultimate) == 100
            @test MortalityTables.ω(cso1980) == 100
        end

        @testset "2001 VBT" begin
            vbt2001 =
                MortalityTables.table("2001 VBT Residual Standard Select and Ultimate - Male Nonsmoker, ANB")


            @test vbt2001.select[35][35] ≈ 0.00036
            @test vbt2001.select[35][35 + 60] ≈ 0.24298
            @test vbt2001.select[95][95] ≈ 0.23815
            @test vbt2001.ultimate[95] ≈ 0.24298
            @test vbt2001.ultimate[120] ≈ 1.0
            @test vbt2001.select[85][120] ≈ 1.0
            @test_throws BoundsError vbt2001.select[150]
            @test omega(vbt2001.select[20]) == 120
            @test omega(vbt2001.ultimate) == 120

            # this tests when there are ultimate rates available
            # but no select rates for the issue age
            @test_throws BoundsError vbt2001.select[116][116] ≈ 0.79988

        end

        @testset "2001 CSO" begin
            cso2001 =
                MortalityTables.table("2001 CSO Super Preferred Select and Ultimate - Male Nonsmoker, ANB")

            @test cso2001.select[35][35] ≈ 0.00037
            @test cso2001.select[35][35 + 60] ≈ 0.26719
            @test cso2001.ultimate[16] ≈ 0.00041
            @test cso2001.ultimate[95] ≈ 0.26719
            @test cso2001.ultimate[120] ≈ 1.0
            @test cso2001.select[85][120] ≈ 1.0
            @test ismissing(cso2001.select[15][15]) # age before table defines rates
            @test_throws BoundsError cso2001.select[150]
            @test_throws BoundsError cso2001.select[70][150]
            @test omega(cso2001.select[16]) == 120
            @test omega(cso2001.select[98]) == 120 # issue #71
            @test omega(cso2001.ultimate) == 120

            @testset "duration ranges" begin
                # select only
                @test all(cso2001.select[29][29:31] .== [0.00029, 0.00035, 0.0004])
                @test all(cso2001.select[35][35:59] .== [
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
                ])

                # crosses into ultimate
                @test all(cso2001.select[35][55:64] .== [
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
                ])
                @test all(cso2001.select[35][35:64] .== [
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
                ])
            end

        end

        @testset "value equality" begin
            ult_table_name = "1980 CSO Basic Table – Male, ANB"
            cso1980_a = MortalityTables.table(ult_table_name)
            cso1980_b = MortalityTables.table(ult_table_name)
            @test cso1980_a == cso1980_b

            select_table_name = "2001 CSO Super Preferred Select and Ultimate - Male Nonsmoker, ANB"
            cso2001_a = MortalityTables.table(select_table_name)
            cso2001_b = MortalityTables.table(select_table_name)
            @test cso2001_a == cso2001_b
        end

    end
end
