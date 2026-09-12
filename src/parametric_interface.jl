abstract type ParametricMortality end

### Generic Functions

"""
    hazard(model,age)

The force of mortality at `age`. More precisely: the ratio of the probability of failure/death to the survival function.

Every parametric law implements `hazard`; `cumhazard`, `survival`, and `decrement` follow from it.
"""
function hazard end

"""
    cumhazard(model,age)
    cumhazard(model,from_age,to_age)

The cumulative force of mortality at `age`. More precisely: the ratio of the cumulative probability of failure/death to the survival function.

The generic method integrates `hazard` numerically from age zero. A law may override the one-argument form with a closed-form expression; the two-age form is the difference of the one-argument form, so a closed form is used automatically wherever it exists.
"""
function cumhazard end

"""
    μ(;m::ParametricMortality,age)

``\\mu_x``: Return the force of mortality at the given age.
"""
function μ(m::ParametricMortality, age)
    return hazard(m,age)
end

cumhazard(m::ParametricMortality, to)       = quadgk(age -> hazard(m, age), 0, to)[1]
cumhazard(m::ParametricMortality, from, to) = cumhazard(m, to) - cumhazard(m, from)
survival(m::ParametricMortality, to)        = exp(-cumhazard(m, to))
survival(m::ParametricMortality, from, to)  = exp(-cumhazard(m, from, to))

# Continuous models need no fractional-age assumption; accept and ignore one.
survival(m::ParametricMortality, to, ::DeathDistribution) = survival(m, to)
survival(m::ParametricMortality, from, to, ::DeathDistribution) = survival(m, from, to)

# A parametric model has no last defined age.
omega(::ParametricMortality) = Inf

(m::ParametricMortality)(x) = μ(m, x)
Base.getindex(m::ParametricMortality,x) = m(x)
Base.broadcastable(pm::ParametricMortality) = Ref(pm)
