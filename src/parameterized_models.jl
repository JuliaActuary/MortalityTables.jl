abstract type ParametricMortality end
"""
    Makeham(a,b,c)

Construct a mortality model following Makeham's law.

"""
struct Makeham <: ParametricMortality
    a
    b
    c
end

"""
    μ(m::ParametricMortality,age)

``\\mu_x``: Return the force of mortality at the given age. 
"""
hazard(m::Makeham,age) = m.a*exp(m.b*age) + m.c
cumhazard(m::Makeham,age) = m.a / m.b * (exp(m.b*age) - 1) + age * m.c
survivorship(m::Makeham,age) = exp(-cumhazard(m,age))
survivorship(m::Makeham,from,to) = survivorship(m::Makeham,to) / survivorship(m::Makeham,from)



function μ(m::ParametricMortality, age) 
    return hazard(m,age)
end

"""
GompertzMakeham(a,b,c)

Construct a mortality model following the full version of GompertzMakeham's law. 
"""
struct GompertzMakeham <: ParametricMortality
    a
    b
    c
end

function μ(m::GompertzMakeham, age) 
    (a * exp(b *age) + c) * exp(-c*x - a/b * (exp(b*x - 1)))
end

"""
Gompertz(b,c)

Construct a mortality model following Gompertz' law of mortality.

This is a special case of Makeham's law.

Calling this will create a `Makeham` model.

"""
function Gompertz(a, b) 
    return Makeham(a, b, 0)
end

# use the integral to calculate the one-year survival
function survivorship(m::ParametricMortality, from_age, to_age) 
    if from_age == to_age
        return 1.0
    else
        return exp(-quadgk(age->μ(m, age), from_age, to_age)[1])
    end
end
survivorship(m::ParametricMortality,to_age) = survivorship(m, 0, to_age)

decrement(m::ParametricMortality,from_age,to_age) = 1 - survivorship(m, from_age, to_age)
decrement(m::ParametricMortality,to_age) = 1 - survivorship(m, to_age)

(m::ParametricMortality)(x) = μ(m, x)
Base.getindex(m::ParametricMortality,x) = m(x)
Base.broadcastable(pm::ParametricMortality) = Ref(pm)