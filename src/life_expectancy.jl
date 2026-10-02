"""
    life_expectancy(table,age)
    life_expectancy(table,age,DeathDistribution)

Calcuate the remaining life expectancy. Assumes curtate life expectancy for tables if not Parametric or DeathDistribution given.

The life_expectancy of the last age defined in the table is set to be `0.0`, even if the table does not end with a rate of `1.0`. An age outside the table is a `BoundsError`.

For a parametric model, survival is integrated from `age` to the model's `omega`: `Inf` for most laws, but finite for a law whose formula ends, such as `Wittstein`. This is the expected remaining lifetime capped at `omega`. If survival is still positive just before `omega`, the remaining probability of death is placed at `omega`, the same convention as a table's last age, whose life_expectancy is `0` even if its last rate is below `1`.

Parametric models accept a `DeathDistribution` and ignore it, since they are continuous.
"""
function life_expectancy(table,age)
    # curtate: sum of the survival probabilities to each later whole age, accumulated in a
    # single pass rather than recomputed from `age` each time. `table[age]` is read first, so
    # an age outside the table is a BoundsError rather than an empty sum. The sum starts from the
    # rates' numeric type, not from `p`, whose rate can be `missing` at the last age.
    p = one(_rate_type(table)) * (1 - table[age])
    s = zero(_rate_type(table))
    for a in (age + 1):lastindex(table)
        s += p
        p *= 1 - table[a]
    end
    return s
end

function life_expectancy(table,age,dist)
    if age == lastindex(table)
        # the type of the integral below, whose nodes are floating-point ages
        return zero(_survival_type(table, age, float(age)))
    else
        QuadGK.quadgk(to -> survival(table,age,to+age,dist),0,lastindex(table)-age)[1]
    end

end

function life_expectancy(table::ParametricMortality,age)
    # over the remaining lifetime t, as for an unbounded law before omega was used; integrating over
    # attained ages [age, omega] gives the same bits but is about 10% slower
    return QuadGK.quadgk(t -> survival(table, age, age + t), 0, omega(table) - age)[1]
end

# continuous models need no fractional-age assumption; accept and ignore one
life_expectancy(table::ParametricMortality, age, ::DeathDistribution) = life_expectancy(table, age)
