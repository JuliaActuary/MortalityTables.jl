abstract type ParametricMortality end

### Generic Functions

"""
    μ(;m::ParametricMortality,age)

``\\mu_x``: Return the force of mortality at the given age. 
"""
function μ(m::ParametricMortality, age) 
    return hazard(m,age)
end

survival(m::ParametricMortality,to_age) = exp(-quadgk(age->μ(m, age), 0, to_age)[1])
survival(m::ParametricMortality,from,to) = survival(m,to) / survival(m,from)

# Continuous models need no fractional-age assumption; accept and ignore one.
survival(m::ParametricMortality, to, ::DeathDistribution) = survival(m, to)
survival(m::ParametricMortality, from, to, ::DeathDistribution) = survival(m, from, to)

# A parametric model has no last defined age.
omega(::ParametricMortality) = Inf

(m::ParametricMortality)(x) = μ(m, x)
Base.getindex(m::ParametricMortality,x) = m(x)
Base.broadcastable(pm::ParametricMortality) = Ref(pm)
