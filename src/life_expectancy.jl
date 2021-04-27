"""
    life_expectancy(table,age)
    life_expectancy(table,age,DeathDistribution)

Calcuate the remaining life expectancy. Assumes curtate life expectancy for tables if not Parametric or DeathDistribution given.
"""
function life_expectancy(table,age)
    sum(survival(table,age,age + dur) for dur in 1:lastindex(table) -age)
end

function life_expectancy(table,age,dist)
    QuadGK.quadgk(to -> survival(table,age,to+age,dist),0,lastindex(table)-age)[1]
end

function life_expectancy(table::ParametricMortality,age)
    QuadGK.quadgk(to -> survival(table,age,to+age),0,Inf)[1]
end