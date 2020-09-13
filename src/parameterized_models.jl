# Reference/Adapted from:
# https://www.mortality.org/Public/HMD_4th_Symposium/Pascariu_poster.pdf
# https://github.com/mpascariu/MortalityLaws/blob/master/R/MortalityLaw_models.R

abstract type ParametricMortality end
"""
    Makeham(;a,b,c)

Construct a mortality model following Makeham's law.

Default args:
    
    a = 0.0002
    b = 0.13
    c = 0.001

"""
Base.@kwdef struct Makeham <: ParametricMortality
    a = 0.0002
    b = 0.13
    c = 0.001
end


function hazard(m::Makeham,age) 
    @unpack a,b,c = m
    return a*exp(b*age) + c
end

function cumhazard(m::Makeham,age) 
    @unpack a,b,c = m
    return a / b * (exp(b*age) - 1) + age * c
end

survivorship(m::Makeham,age) = exp(-cumhazard(m,age))


"""
    Gompertz(a,b)

Construct a mortality model following Gompertz' law of mortality.

This is a special case of Makeham's law and will `Makeham` model where `c=0`.

Default args:

    a = 0.0002
    b = 0.13

"""
function Gompertz(;a=0.0002, b=0.13) 
    return Makeham(a=a, b=b, c=0)
end

"""
    InverseGompertz(;a,b,c)

Construct a mortality model following InverseGompertz's law.

Default args:
    
    m = 49
    σ = 7.7

"""
Base.@kwdef struct InverseGompertz <: ParametricMortality
    m = 49
    σ = 7.7
end


function hazard(m::InverseGompertz,age)
    @unpack m,σ = m 
    return 1 / σ * exp(-(age - m)/σ) / (exp(exp(-(age - m)/σ)) - 1)
end

cumhazard(m::InverseGompertz,age) = -log(survivorship(m,age))

function survivorship(m::InverseGompertz,age) 
    @unpack m,σ = m 
    return (1 - exp(-exp(-(age - m)/σ))) / (1 - exp(-exp(m/σ)))
end

"""
    Opperman(a,b,c)

Construct a mortality model following Opperman's law of mortality.

Default args:

    a = 0.04
    b = 0.0004
    c = 0.001
"""
Base.@kwdef struct Opperman <: ParametricMortality
    a = 0.04
    b = 0.0004
    c = 0.001
end

function hazard(m::Opperman,age) 
    @unpack a,b,c = m
    return max(a / √(age+1) - b + c * √(age+1),0.0)
end

"""
    Thiele(a,b,c,d,e,f,g)

Construct a mortality model following Opperman's law of mortality.

Default args:

    a = 0.02474 
    b = 0.3
    c = 0.004
    d = 0.5
    e = 25
    f = 0.0001
    g = 0.13
"""
Base.@kwdef struct Thiele <: ParametricMortality
    a = 0.02474 
    b = 0.3
    c = 0.004
    d = 0.5
    e = 25
    f = 0.0001
    g = 0.13
end

function hazard(m::Thiele,age) 
    @unpack a,b,c,d,e,f,g = m
    μ₁ = a * exp(-b * age)
    μ₂ = c * exp(-0.5 * d * (age - e)^2)
    μ₃ = f * exp(g * age)

    if age == 0 
        return μ₁ + μ₃
    else
        return  μ₁ + μ₂ + μ₃
    end
end

"""
    Wittstein(a,b,m,n)

Construct a mortality model following Wittstein's law of mortality.

Default args:

    a = 1.5
    b = 1.
    n = 0.5
    m = 100

"""
Base.@kwdef struct Wittstein <: ParametricMortality
    a = 1.5
    b = 1.
    n = 0.5
    m = 100
end

function hazard(m::Wittstein,age)
    @unpack a,b,m,n = m
    return (1/b) * a ^ -((b * age) ^ n) + a^ -((m -  age) ^ n) 
end

"""
    Weibull(m,σ)

Construct a mortality model following Weibull's law of mortality.

Note that if σ > m, then the mode of the density is 0 and hx is a non-increasing function of x, while if σ < m, then the mode is greater than 0 and hx is an increasing function. 
 - `m >0` is a measure of location
 - `σ >0` is measure of dispersion


 Default args:

    m = 1
    σ = 2
"""
Base.@kwdef struct Weibull <: ParametricMortality
    m = 1
    σ = 2
end

function hazard(m::Weibull,age)
    @unpack m,σ =m
    if age == 0
        return 1.0
    else 
        return 1 / σ * (age / m)^(m / σ - 1)
    end
end

function cumhazard(m::Weibull,age)
    @unpack m,σ =m
    return (age / m) ^ (m / σ)
end

function survivorship(m::Weibull,age) 
    return exp(-cumhazard(m,age))
end

"""
    InverseWeibull(m,σ)

Construct a mortality model following Weibull's law of mortality.

The Inverse-Weibull proves useful for modelling the childhood and teenage years, because the logarithm of h(x) is a concave function.
 - `m >0` is a measure of location
 - `σ >0` is measure of dispersion

 Default args:

    m = 5
    σ = 10

"""
Base.@kwdef struct InverseWeibull <: ParametricMortality
    m = 5
    σ = 10
end

function hazard(m::InverseWeibull,age)
    @unpack m,σ = m
    return (1/σ) * (age/m)^(-m/σ - 1) / (exp((age/m)^(-m/σ)) - 1)
end

function cumhazard(m::InverseWeibull,age)
    @unpack m,σ = m
    return -log(1 - exp(-(age/m)^(-m/σ)))
end

function survivorship(m::InverseWeibull,age) 
    return exp(-cumhazard(m,age))
end

"""
    Perks(a,b,c,d)

Construct a mortality model following Perks' law of mortality.

Default args:

    a = 0.002
    b = 0.13
    c = 0.01
    d = 0.01
"""
Base.@kwdef struct Perks <: ParametricMortality
    a = 0.002
    b = 0.13
    c = 0.01
    d = 0.01
end

function hazard(m::Perks,age) 
    @unpack a,b,c,d = m
    return (a + b*c^age) / (b*(c^-age) + 1 + d*c^age)
end

"""
    VanderMaen(a,b,c,i,n)

Construct a mortality model following VanderMaen's law of mortality.

Default args:
    
    a = 0.01
    b = 1
    c = 0.01
    i = 100
    n = 200
"""
Base.@kwdef struct VanderMaen <: ParametricMortality
    a = 0.01
    b = 1.
    c = 0.01
    i = 100.
    n = 200.
end

function hazard(m::VanderMaen,age)
    @unpack a,b,c,i,n = m
    return a + b*age + c*(age^2) + i/(n - age)
end

"""
    VanderMaen2(a,b,i,n)

Construct a mortality model following VanderMaen2's law of mortality.

Default args:

    a = 0.01
    b = 1
    i = 100
    n = 200

"""
Base.@kwdef struct VanderMaen2 <: ParametricMortality
    a = 0.01
    b = 1.
    i = 100.
    n = 200.
end

function hazard(m::VanderMaen2,age)
    @unpack a,b,i,n = m 
    return a + b * age + i/(n - age)
end

"""
    StrehlerMildvan(k,v₀,b,d)

Construct a mortality model following StrehlerMildvan's law of mortality.

Default args:

    k   = 0.01
    v₀  = 2.5
    b   = 0.2
    d   = 6.0

"""
Base.@kwdef struct StrehlerMildvan <: ParametricMortality
    k   = 0.01
    v₀  = 2.5
    b   = 0.2
    d   = 6.0
end

function hazard(m::StrehlerMildvan,age)
    @unpack k,v₀,b,d = m
    return  k * exp(-v₀ * (1 - b * age) / d)
end

"""
    Beard(a,b,k)

Construct a mortality model following Beard's law of mortality.

Default args:
    
    a = 0.002
    b = 0.13
    k = 1.
"""
Base.@kwdef struct Beard <: ParametricMortality
    a = 0.002
    b = 0.13
    k = 1.
end

function hazard(m::Beard,age)
    @unpack a,b,k = m
    return  a * exp(b*age) / (1 + k * a * exp(b*age))
end

"""
    MakehamBeard(a,b,c,k)

Construct a mortality model following MakehamBeard's law of mortality.

Default args:

    a = 0.002
    b = 0.13
    c = 0.01
    k = 1.
"""
Base.@kwdef struct MakehamBeard <: ParametricMortality
    a = 0.002
    b = 0.13
    c = 0.01
    k = 1.
end

function hazard(m::MakehamBeard,age)
    @unpack a,b,c,k = m
    return  a * exp(b*age) / (1 + k * a * exp(b*age)) + c
end

"""
    Quadratic(a,b,c)

Construct a mortality model following Quadratic law of mortality.

Default args:

    a = 0.01
    b = 1.
    c = 0.01
"""
Base.@kwdef struct Quadratic <: ParametricMortality
    a = 0.01
    b = 1.
    c = 0.01
end

function hazard(m::Quadratic,age)
    @unpack a,b,c = m
    return  a + b * age + c * age^2
end

"""
    GammaGompertz(a,b,γ)

Construct a mortality model following GammaGompertz law of mortality.

Default args:

    a = 0.002
    b = 0.13
    γ = 1
"""
Base.@kwdef struct GammaGompertz <: ParametricMortality
    a = 0.002
    b = 0.13
    γ = 1
end

function hazard(m::GammaGompertz,age)
    @unpack a,b,γ = m
    return  (a * exp(b * age)) / (1 + ( a * γ / b) * (exp(b * age) - 1))
end

"""
    Siler(a,b,c,d,e)

Construct a mortality model following Siler law of mortality.

Default args:

    a = 0.0002
    b = 0.13
    c = 0.001
    d = 0.001
    e = 0.013
"""
Base.@kwdef struct Siler <: ParametricMortality
    a = 0.0002
    b = 0.13
    c = 0.001
    d = 0.001
    e = 0.013
end

function hazard(m::Siler,age)
    @unpack a,b,c,d,e = m
    return  a * exp(-b* age) + c + d * exp(e * age)
end

"""
    HeligmanPollard(a,b,c,d,e,f,g,h)

Construct a mortality model following HeligmanPollard law of mortality with 8 parameters.

Default args:

    a = 0.0002
    b = 0.13
    c = 0.001
    d = 0.001
    e = 0.013
"""
Base.@kwdef struct HeligmanPollard <: ParametricMortality
    a = .0005
    b = .004
    c = .08
    d = .001
    e = 10
    f = 17
    g = .00005
    h = 1.1
end

function hazard(m::HeligmanPollard,age)
    @unpack a,b,c,d,e,f,g,h = m
    μ₁ = a^((age + b)^c) + g * h^age
    μ₂ = d * exp(-e * (log(age/f))^2)
    η = age == 0 ?  μ₁ :  μ₁ + μ₂
    return  η / (1 + η)
end

"""
    HeligmanPollard2(a,b,c,d,e,f,g,h)

Construct a mortality model following HeligmanPollard (alternate) law of mortality with 8 parameters.

Default args:

    a = .0005
    b = .004
    c = .08
    d = .001
    e = 10
    f = 17
    g = .00005
    h = 1.1
"""
Base.@kwdef struct HeligmanPollard2 <: ParametricMortality
    a = 0.0005
    b = 0.004
    c = 0.08
    d = 0.001
    e = 10.
    f = 17.
    g = 0.00005
    h = 1.1
end

function hazard(m::HeligmanPollard2,age)
    @unpack a,b,c,d,e,f,g,h = m
    μ₁ = a^((age + b)^c) + (g * h^age) / (1 + g * h ^ age)
    μ₂ = d * exp(-e * (log(age/f))^2)
    return age == 0 ?  μ₁ :  μ₁ + μ₂
end

"""
    HeligmanPollard3(a,b,c,d,e,f,g,h,k)

Construct a mortality model following HeligmanPollard (alternate) law of mortality with 9 parameters.

Default args:

    a = .0005
    b = .004
    c = .08
    d = .001
    e = 10
    f = 17
    g = .00005
    h = 1.1
    k= 1.
"""
Base.@kwdef struct HeligmanPollard3 <: ParametricMortality
    a = .0005
    b = .004
    c = .08
    d = .001
    e = 10
    f = 17
    g = .00005
    h = 1.1
    k = 1.
end

function hazard(m::HeligmanPollard3,age)
    @unpack a,b,c,d,e,f,g,h,k = m
    μ₁ = a^((age + b)^c) + (g * h^age) / (1 + k * g * h ^ age)
    μ₂ = d * exp(-e * (log(age/f))^2)
    return age == 0 ?  μ₁ :  μ₁ + μ₂
end

"""
    HeligmanPollard4(a,b,c,d,e,f,g,h,k)

Construct a mortality model following HeligmanPollard (alternate) law of mortality with 9 parameters.

Default args:

    a = .0005
    b = .004
    c = .08
    d = .001
    e = 10
    f = 17
    g = .00005
    h = 1.1
    k= 1.
"""
Base.@kwdef struct HeligmanPollard4 <: ParametricMortality
    a = .0005
    b = .004
    c = .08
    d = .001
    e = 10
    f = 17
    g = .00005
    h = 1.1
    k = 1.
end

function hazard(m::HeligmanPollard4,age)
    @unpack a,b,c,d,e,f,g,h,k = m
    μ₁ = a^((age + b)^c) + (g * h^(age ^ k)) / (1 + g * h ^ (age ^ k))
    μ₂ = d * exp(-e * (log(age/f))^2)
    return age == 0 ?  μ₁ :  μ₁ + μ₂
end

"""
    RogersPlanck(a₀, a₁, a₂, a₃, a, b, c, d, u)

Construct a mortality model following RogersPlanck law of mortality.

Default args:

    a₀ = 0.0001
    a₁ = 0.02
    a₂ = 0.001
    a₃ = 0.0001
    a  = 2.
    b  = 0.001
    c  = 100.
    d  = 0.1
    u  = 0.33

"""
Base.@kwdef struct RogersPlanck <: ParametricMortality
    a₀ = 0.0001
    a₁ = 0.02
    a₂ = 0.001
    a₃ = 0.0001
    a  = 2.
    b  = 0.001
    c  = 100.
    d  = 0.1
    u  = 0.33
end

function hazard(m::RogersPlanck,age) 
    @unpack a₀,a₁,a₂,a₃,a,b,c,d,u = m
    return  a₀ + a₁ * exp(-a * age) + a₂ * exp(b*(age - u) - exp(-c*(age - u))) + a₃*exp(d*age)
end


"""
    Martinelle(a,b,c,d,k)

Construct a mortality model following Martinelle's law of mortality.

Default args:

    a = 0.001
    b = 0.13
    c = 0.001
    d = 0.1
    k = 0.001
"""
Base.@kwdef struct Martinelle <: ParametricMortality
    a = 0.001
    b = 0.13
    c = 0.001
    d = 0.1
    k = 0.001
end

function hazard(m::Martinelle,age)
    @unpack a,b,c,d,k = m
    return  (a*exp(b*age) + c) / (1 + d*exp(b * age)) + k*exp(b * age)
end


"""
    Kostaki(a,b,c,d,e1,e2,f,g,h)

Construct a mortality model following Kostaki's law of mortality. A nine-parameter adaptation of `HeligmanPollard`.

 > Kostaki, A. (1992). A nine‐parameter version of the Heligman‐Pollard formula. Mathematical Population Studies, 3(4), 277–288. doi:10.1080/08898489209525346 

Default args:

    a = 0.0005
    b = 0.01
    c = 0.10
    d = 0.001
    e1 = 3.
    e2 = 0.1
    f = 25.
    g = .00005
    h = 1.1


"""
Base.@kwdef struct Kostaki <: ParametricMortality
    a = 0.0005
    b = 0.01
    c = 0.10
    d = 0.001
    e1 = 3.
    e2 = 0.1
    f = 25.
    g = .00005
    h = 1.1
end

function hazard(m::Kostaki,age) 
    @unpack a,b,c,d,e1,e2,f,g,h = m
    μ₁ = a^((age + b)^c) + g*h^age 
    if age <= f
        μ₂ =d *exp(-(e1*log(age/f))^2)
    else
        μ₂ =d *exp(-(e2*log(age/f))^2)
    end
    
    η = age == 0 ? μ₁ : μ₁ + μ₂

    return η / (1+η)
end

"""
    Kannisto(a,b)

Construct a mortality model following Kannisto's law of mortality.

Default args:

    a = 0.5
    b = 0.13
"""
Base.@kwdef struct Kannisto <: ParametricMortality
    a = 0.5
    b = 0.13
end

function hazard(m::Kannisto,age)
    @unpack a,b = m
    return  a * exp(b * age) / (1 + a * exp(b*age))
end

function cumhazard(m::Kannisto,age)
    @unpack a,b = m
    return  1/a * log((1 + b*exp(b*age)) / (1 + a))
end

function  survivorship(m::Kannisto,age)
    return exp(-cumhazard(m,age))
end


"""
    KannistoMakeham(a,b)

Construct a mortality model following KannistoMakeham's law of mortality.

Default args:

    a = 0.5
    b = 0.13
    c = 0.001
"""
Base.@kwdef struct KannistoMakeham <: ParametricMortality
    a = 0.5
    b = 0.13
    c = 0.001
end

function hazard(m::KannistoMakeham,age)
    @unpack a,b,c = m
    return  a * exp(b * age) / (1 + a * exp(b*age)) + c
end

### Generic Functions

"""
    μ(m::ParametricMortality,age)

``\\mu_x``: Return the force of mortality at the given age. 
"""
function μ(m::ParametricMortality, age) 
    return hazard(m,age)
end


# # use the integral to calculate the one-year survival
# function survivorship(m::ParametricMortality, from_age, to_age) 
#     if from_age == to_age
#         return 1.0
#     else
#         return exp(-quadgk(age->μ(m, age), from_age, to_age)[1])
#     end
# end
# survivorship(m::ParametricMortality,to_age) = survivorship(m, 0, to_age)

survivorship(m::ParametricMortality,to_age) = exp(-quadgk(age->μ(m, age), 0, to_age)[1])
survivorship(m,from,to) = survivorship(m,to) / survivorship(m,from)
decrement(m::ParametricMortality,from_age,to_age) = 1 - survivorship(m, from_age, to_age)
decrement(m::ParametricMortality,to_age) = 1 - survivorship(m, to_age)

(m::ParametricMortality)(x) = μ(m, x)
Base.getindex(m::ParametricMortality,x) = m(x)
Base.broadcastable(pm::ParametricMortality) = Ref(pm)