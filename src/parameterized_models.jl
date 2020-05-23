
"""
    Makeham(a,b,c)

Construct a mortality model following Makeham's law:

``\\mu_x = A + Bc^x``

"""
struct Makeham
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
μ(m::Makeham,age) = m.a + m.b * m.c ^ age

# use the integral to calculate the one-year survival
"""
    p(m::Makeham,age,time) 
    p(m::Makeham,age) 

With `time` argument: ``_tp_x``, the survival from `age` to `age + time`.
Otherwise: ``p_x``, the survival from `age` to `age  + 1`.
"""
p(m::Makeham,x,t) = exp(-quadgk(s->μ(m,x+s),0,t)[1])
p(m::Makeham,x) = p(m::Makeham,x,1)

"""
    q(m::Makeham,age,time) 
    q(m::Makeham,age) 

With `time` argument: ``_tq_x``, the cumulative mortality from `age` to `age + time`.
Otherwise: ``q_x``, the cumulative mortality from `age` to `age  + 1`.
"""
q(m::Makeham,x,t) = 1 - p(m,x,t)
q(m::Makeham,x) = 1 - p(m,x,1)

