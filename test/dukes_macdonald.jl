using MortalityTables: SelectUltimateTable, TableMetaData, table, dukes_macdonald1, dukes_macdonald2, dukes_macdonald3

@testset "Dukes/MacDonald" begin
  base_lapses = 0.10
  total_lapses = 0.85
  effectiveness = 0.80
  @testset "from newsletter Doll, 2003" begin
    select = 0.01
    point_in_scale = 0.03
    issue_age = 1
    attained_age = 2
    testSelectUltimate = SelectUltimateTable(
      [[0.005, point_in_scale], [nothing, select]], # nothing is like offset array
      [0.02, 0.04, 0.05],
      TableMetaData()
    )
    @test 3.67 == round(dukes_macdonald1(testSelectUltimate, issue_age, attained_age, base_lapses, total_lapses, effectiveness) / point_in_scale; digits=2)
    @test 2.33 == round(dukes_macdonald2(testSelectUltimate, issue_age, attained_age, base_lapses, total_lapses, effectiveness) / point_in_scale; digits=2)
    @test 2.00 == round(dukes_macdonald3(testSelectUltimate, issue_age, attained_age, base_lapses, total_lapses, effectiveness) / point_in_scale; digits=2)
  end
  @testset "using real data" begin
    issue_age = 50
    attained_age = 60
    tbl = table(3299)
    dukes_methods = [
      dukes_macdonald1,
      dukes_macdonald2,
      dukes_macdonald3
    ]
    results = [
      dukes_method(tbl, issue_age, attained_age, base_lapses, total_lapses, effectiveness) / tbl.select[attained_age][attained_age]
      for dukes_method in dukes_methods
    ]
    @test 1 < results[3] < results[2] < results[1] < 8
  end
end