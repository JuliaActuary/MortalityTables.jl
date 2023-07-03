using MortalityTables: SelectUltimateTable

"""
  dukes_macdonald1(point_in_scale, select_rate, base_lapses, total_lapses, effectiveness)

  dukes_macdonald1(selectultimate, issue_age, attained_age, base_lapses, total_lapses, effectiveness)

Two methods are available, one takes the point_in_scale and select_rate directly. See https://www.soa.org/globalassets/assets/library/newsletters/product-development-news/2003/july/pdn-2003-iss56-doll-a.pdf for information on terminology. 
The other method takes a SelectUltimateTable from MortalityTable.jl and calculates the point_in_scale and select_rate from the table.
"""
function dukes_macdonald1(point_in_scale, select_rate, base_lapses, total_lapses, effectiveness)
  excess_lapses = total_lapses - base_lapses
  select_excess_lapses = effectiveness * excess_lapses
  not_lapses = 1 - total_lapses
  deteriorated_rate = (point_in_scale * (select_excess_lapses + not_lapses) - select_excess_lapses * select_rate) / not_lapses
  return deteriorated_rate
end

function dukes_macdonald1(selectultimate, issue_age, attained_age, base_lapses, total_lapses, effectiveness)
  point_in_scale = selectultimate.select[issue_age][attained_age]
  select_rate = selectultimate.select[attained_age][attained_age]
  return dukes_macdonald1(point_in_scale, select_rate, base_lapses, total_lapses, effectiveness)
end

"""
  dukes_macdonald2(point_in_scale, select_rate, base_lapses, total_lapses, effectiveness)
  
  dukes_macdonald2(selectultimate, issue_age, attained_age, base_lapses, total_lapses, effectiveness)

Two methods are available, one takes the point_in_scale and select_rate directly. See https://www.soa.org/globalassets/assets/library/newsletters/product-development-news/2003/july/pdn-2003-iss56-doll-a.pdf for information on terminology. 
The other method takes a SelectUltimateTable from MortalityTable.jl and calculates the point_in_scale and select_rate from the table.
"""
function dukes_macdonald2(point_in_scale, select_rate, base_lapses, total_lapses, effectiveness)
  excess_lapses = total_lapses - base_lapses
  select_excess_lapses = effectiveness * excess_lapses
  nonselect_excess_lapses = excess_lapses - select_excess_lapses
  not_lapses = 1 - total_lapses
  deteriorated_rate = (point_in_scale * (not_lapses + excess_lapses) - select_excess_lapses * select_rate) / (nonselect_excess_lapses + not_lapses)
  return deteriorated_rate
end

function dukes_macdonald2(selectultimate, issue_age, attained_age, base_lapses, total_lapses, effectiveness)
  point_in_scale = selectultimate.select[issue_age][attained_age]
  select_rate = selectultimate.select[attained_age][attained_age]
  return dukes_macdonald2(point_in_scale, select_rate, base_lapses, total_lapses, effectiveness)
end

"""
  dukes_macdonald3(point_in_scale, select_rate, base_lapses, total_lapses, effectiveness)
  
  dukes_macdonald3(selectultimate, issue_age, attained_age, base_lapses, total_lapses, effectiveness)

Two methods are available, one takes the point_in_scale and select_rate directly. See https://www.soa.org/globalassets/assets/library/newsletters/product-development-news/2003/july/pdn-2003-iss56-doll-a.pdf for information on terminology. 
The other method takes a SelectUltimateTable from MortalityTable.jl and calculates the point_in_scale and select_rate from the table.
"""
function dukes_macdonald3(point_in_scale, select_rate, base_lapses, total_lapses, effectiveness)
  excess_lapses = total_lapses - base_lapses
  select_excess_lapses = effectiveness * excess_lapses
  deteriorated_rate = (point_in_scale - select_excess_lapses * select_rate) / (1 - select_excess_lapses)
  return deteriorated_rate
end

function dukes_macdonald3(selectultimate, issue_age, attained_age, base_lapses, total_lapses, effectiveness)
  point_in_scale = selectultimate.select[issue_age][attained_age]
  select_rate = selectultimate.select[attained_age][attained_age]
  return dukes_macdonald3(point_in_scale, select_rate, base_lapses, total_lapses, effectiveness)
end