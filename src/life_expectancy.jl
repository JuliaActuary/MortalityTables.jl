"""
    life_expectancy(table,age)
    life_expectancy(table,age,DeathDistribution)

Calcuate the remaining life expectancy. Assumes curtate life expectancy for tables if not Parametric or DeathDistribution given.

The life_expectancy of the last age defined in the table is set to be `0.0`, even if the table does not end with a rate of `1.0`.
"""
function life_expectancy(table,age)
    if age == lastindex(table)
        return 0.
    else
        sum(survival(table,age,age + dur) for dur in 1:lastindex(table) -age)
    end
end

function life_expectancy(table,age,dist)
    if age == lastindex(table)
        return 0.
    else
        QuadGK.quadgk(to -> survival(table,age,to+age,dist),0,lastindex(table)-age)[1]
    end

end

function life_expectancy(table::ParametricMortality,age)
    QuadGK.quadgk(to -> survival(table,age,to+age),0,Inf)[1]
end

life_expectancy(table::MortalityTable,args...) = throw(ArgumentError("The first argument should be a vector of rates instead of an entire table. E.g. `table.ulitmate` or `table.select[age]`."))