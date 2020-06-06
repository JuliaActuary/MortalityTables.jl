abstract type ParametricMortality end
"""
    Makeham(a,b,c)

Construct a mortality model following Makeham's law:

``\\mu_x = A + Bc^x``

"""
struct Makeham <: ParametricMortality
    a
    b
    c
end

"""
    Gompertz(b,c)

Construct a mortality model following Gompertz' law of mortality:

``\\mu_x = Bc^x``

This is a special case of Makeham's law ``\\mu_x = A + Bc^x``, where ``A=0``.

Calling this will create a `Makeham` model that can be called with `q` or `p`.

"""
function Gompertz(b,c) 
    return Makeham(0,b,c)
end

"""
    μ(m::Makeham,age)

``\\mu_x``: Return the force of mortality at the given age. 
"""
μ(m::Makeham,age::UnitRange) = m.a .+ m.b .* m.c .^ collect(age)
function μ(m::Makeham,age) 
    m.a .+ m.b * m.c .^ age
end

# use the integral to calculate the one-year survival
"""
    p(m::Makeham,age,time) 
    p(m::Makeham,age) 

With `time` argument: ``_tp_x``, the survival from `age` to `age + time`.
Otherwise: ``p_x``, the survival from `age` to `age  + 1`.
"""
function survivorship(m::Makeham,from_age,to_age) 
    if from_age == to_age
        return 1.0
    else
        return exp(-quadgk(age->μ(m,age),from_age,to_age)[1])
    end
end
survivorship(m::Makeham,to_age) = survivorship(m::Makeham,0,to_age)

"""
    q(m::Makeham,age,time) 
    q(m::Makeham,age) 

With `time` argument: ``_tq_x``, the cumulative mortality from `age` to `age + time`.
Otherwise: ``q_x``, the cumulative mortality from `age` to `age  + 1`.
"""
cumulative_decrement(m::Makeham,from_age,to_age) = 1 - survivorship(m,from_age,to_age)
cumulative_decrement(m::Makeham,to_age) = 1 - survivorship(m,to_age)

(m::ParametricMortality)(x) = μ(m,x)
Base.getindex(m::ParametricMortality,x) = m(x)