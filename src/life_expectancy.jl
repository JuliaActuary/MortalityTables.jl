"""
    life_expectancy(table,age)
    life_expectancy(table,age,DeathDistribution)

Calcuate the remaining life expectancy. Assumes curtate life expectancy for tables if not Parametric or DeathDistribution given.

The life_expectancy of the last age defined in the table (or any age past it) is set to be `0.0`, even if the table does not end with a rate of `1.0`.

Parametric models accept a `DeathDistribution` and ignore it, since they are continuous.
"""
function life_expectancy(table,age)
    # curtate: sum of the survival probabilities to each later whole age,
    # accumulated in a single pass rather than recomputed from `age` each time
    s = 0.0
    p = 1.0
    for a in age:lastindex(table)-1
        p *= 1 - table[a]
        s += p
    end
    return s
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

# continuous models need no fractional-age assumption; accept and ignore one
life_expectancy(table::ParametricMortality, age, ::DeathDistribution) = life_expectancy(table, age)

life_expectancy(table::MortalityTable,args...) = throw(ArgumentError("The first argument should be a vector of rates instead of an entire table. E.g. `table.ulitmate` or `table.select[age]`."))