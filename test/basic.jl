
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


        @test survival(q4, -1) ≈ 1.0
        @test survival(q4, 4, -1) ≈ 1.0
        @test decrement(q4, -1) ≈ 0.0
        @test decrement(q4, 4, -1) ≈ 0.0
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
