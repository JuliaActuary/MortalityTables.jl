abstract type ParametricMortality end

### Generic Functions

"""
    hazard(model,age)

The force of mortality at `age`: the density of the age at death divided by the probability of surviving to `age`, ``\\mu_x = f(x)/S(x) = -S'(x)/S(x)``.

Every parametric law implements `hazard`; `cumhazard`, `survival`, and `decrement` follow from it.
"""
function hazard end

"""
    cumhazard(model,age)
    cumhazard(model,from_age,to_age)

The integrated hazard up to `age`, ``H(x) = \\int_0^x \\mu_s \\, ds``, equivalently ``-\\log S(x)``, so that `survival(model, age) == exp(-cumhazard(model, age))`. The two-age form integrates from `from_age` to `to_age`: ``-\\log(S(\\text{to})/S(\\text{from}))``.

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

# 1 - exp(-H), without rounding a small decrement to zero
decrement(m::ParametricMortality, to)       = -expm1(-cumhazard(m, to))
decrement(m::ParametricMortality, from, to) = -expm1(-cumhazard(m, from, to))

# Continuous models need no fractional-age assumption; accept and ignore one.
survival(m::ParametricMortality, to, ::DeathDistribution) = survival(m, to)
survival(m::ParametricMortality, from, to, ::DeathDistribution) = survival(m, from, to)
decrement(m::ParametricMortality, to, ::DeathDistribution) = decrement(m, to)
decrement(m::ParametricMortality, from, to, ::DeathDistribution) = decrement(m, from, to)

# The last age at which a law is defined: every age, unless the law's formula ends (see the
# `omega` docstring and the methods beside `Wittstein`, `VanderMaen` and `VanderMaen2`).
omega(::ParametricMortality) = Inf

(m::ParametricMortality)(x) = μ(m, x)
Base.broadcastable(pm::ParametricMortality) = Ref(pm)
