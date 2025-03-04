# Reference/Adapted from:
# https://mortality.org/File/GetDocument/Public/HMD_4th_Symposium/Pascariu_poster.pdf
# https://github.com/mpascariu/MortalityLaws/blob/master/R/MortalityLaw_models.R

abstract type ParametricMortality end
"""
    Makeham(;a,b,c)

Construct a mortality model following Makeham's law.

``
\\mathrm{hazard} \\left( {\\rm age} \\right) =  ae^{bx} + c
``

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

"""
    hazard(model,age)

The force of mortality at `age`. More precisely: the ratio of the probability of failure/death to the survival function.
"""
function hazard(m::Makeham,age) 
    @unpack a,b,c = m
    return a*exp(b*age) + c
end

"""
    cumhazard(model,age)

The cumulative force of mortality at `age`. More precisely: the ratio of the cumulative probability of failure/death to the survival function.
"""
function cumhazard(m::Makeham,age) 
    @unpack a,b,c = m
    return a / b * (exp(b*age) - 1) + age * c
end

survival(m::Makeham,age) = exp(-cumhazard(m,age))


"""
    Makeham2(;μ,σ,c)
Construct a mortality model following Makeham's law. Alternative formulation.
``
\\mathrm{hazard} \\left( {\\rm age} \\right) =  (1/sigma)e^{(x-mu)/sigma} + c
``
Default args:
    
    μ = 49
    σ = 7.692308
    c = 0.001
"""
Base.@kwdef struct Makeham2 <: ParametricMortality
    μ = 49
    σ = 7.692308
    c = 0.001
end

"""
    hazard(model,age)
The force of mortality at `age`. More precisely: the ratio of the probability of failure/death to the survival function.
"""
function hazard(m::Makeham2,age) 
    @unpack μ, σ, c = m
    return (1.0 / σ) .* exp.((age .- μ) ./ σ) .+ c
end

"""
    cumhazard(model,age)
The cumulative force of mortality at `age`. More precisely: the ratio of the cumulative probability of failure/death to the survival function.
"""
function cumhazard(m::Makeham2,age) 
    @unpack μ, σ, c = m
    return exp(-μ / σ) .* (exp.(age ./ σ) .- 1.0) .+ age .* c
end

survival(m::Makeham2,age) = exp.(-cumhazard(m,age))

"""
    Gompertz(;a,b)

Construct a mortality model following Gompertz' law of mortality.

``
\\mathrm{hazard} \\left( {\\rm age} \\right) =  ae^{bx}
``

This is a special case of Makeham's law and will `Makeham` model where `c=0`.

Default args:

    a = 0.0002
    b = 0.13

"""
function Gompertz(;a=0.0002, b=0.13) 
    return Makeham(a=a, b=b, c=0)
end

"""
    Gompertz2(;μ,σ)
Construct a mortality model following Gompertz' law of mortality. Alternative formulation.
``
\\mathrm{hazard} \\left( {\\rm age} \\right) =  (1/sigma)e^{(x-mu)/sigma}
``
This is a special case of Makeham's law alternative formulation and will `Makeham` model where `c=0`.
Default args:
    μ = 49
    σ = 7.7
"""
Base.@kwdef struct Gompertz2 <: ParametricMortality
    μ = 49
    σ = 7.7
end

"""
    hazard(model,age)
The force of mortality at `age`. More precisely: the ratio of the probability of failure/death to the survival function.
"""
function hazard(m::Gompertz2,age) 
    @unpack μ, σ = m
    return (1.0 / σ) .* exp.((age .- μ) ./ σ)
end

"""
    cumhazard(model,age)
The cumulative force of mortality at `age`. More precisely: the ratio of the cumulative probability of failure/death to the survival function.
"""
function cumhazard(m::Gompertz2,age) 
    @unpack μ, σ = m
    return exp(-μ / σ) .* (exp.(age ./ σ) .- 1.0)
end

survival(m::Gompertz2,age) = exp.(-cumhazard(m,age))

"""
    InverseGompertz(;a,b,c)

Construct a mortality model following InverseGompertz's law.

```math
\\begin{aligned}
\\mathrm{hazard} \\left( {\\rm age} \\right) &= \\frac{1}{\\sigma}e^\\frac{age-m}{\\sigma}/e^{e^\\frac{-(age-m)}{\\sigma}-1}``
\\\\
\\mathrm{survival} \\left( {\\rm age} \\right) &= \\frac{1 - e^{ - e^{\\frac{ - \\left( {\\rm age} - m \\right)}{\\sigma}}}}{1 - e^{ - e^{\\frac{m}{\\sigma}}}}``
\\end{aligned}
```

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

cumhazard(m::InverseGompertz,age) = -log(survival(m,age))

function survival(m::InverseGompertz,age) 
    @unpack m,σ = m 
    return (1 - exp(-exp(-(age - m)/σ))) / (1 - exp(-exp(m/σ)))
end

"""
    Opperman(;a,b,c)

Construct a mortality model following Opperman's law of mortality.

``
\\mathrm{hazard} \\left( {\\rm age} \\right) = \\frac{a}{\\sqrt{age}} + b +c\\sqrt[3]{age}
``

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
    Thiele(;a,b,c,d,e,f,g)

Construct a mortality model following Opperman's law of mortality.

```math
\\begin{aligned}
\\mu_1 &= a \\cdot e^{\\left(  - b \\right) \\cdot {\\rm age}}
\\\\
\\mu_2 &= c \\cdot e^{-0.5 \\cdot d \\cdot \\left( {\\rm age} - e \\right)^{2}}
\\\\
\\mu_3 &= f \\cdot e^{g \\cdot {\\rm age}}
\\\\
\\mathrm{hazard} \\left( {\\rm age} \\right) &= \\begin{cases}
\\mu_1 + \\mu_3 & \\text{if } \\left( {\\rm age} = 0 \\right)\\\\
\\mu_1 + \\mu_2 + \\mu_3 & \\text{otherwise}
\\end{cases}
\\end{aligned}
```
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
    Wittstein(;a,b,m,n)

Construct a mortality model following Wittstein's law of mortality.

``\\mathrm{hazard} \\left( {\\rm age} \\right) = \\frac{1}{b} \\cdot a^{ - \\left( b \\cdot {\\rm age} \\right)^{n}} + a^{ - \\left( m - {\\rm age} \\right)^{n}}``

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
    Weibull(;m,σ)

Construct a mortality model following Weibull's law of mortality.

Note that if σ > m, then the mode of the density is 0 and hx is a non-increasing function of x, while if σ < m, then the mode is greater than 0 and hx is an increasing function. 
 - `m >0` is a measure of location
 - `σ >0` is measure of dispersion

```math
\\begin{aligned}
\\mathrm{hazard} \\left( {\\rm age} \\right) = \\frac{1}{\\sigma} \\cdot \\left( \\frac{{\\rm age}}{m} \\right)^{\\frac{m}{\\sigma} - 1}
\\\\
\\mathrm{cumhazard} \\left( {\\rm age} \\right) = \\left( \\frac{{\\rm age}}{m} \\right)^{\\frac{m}{\\sigma}}
\\\\
\\mathrm{survival} \\left( {\\rm age} \\right) =  e^{ - \\mathrm{cumhazard} \\left( m, {\\rm age} \\right)}
\\end{aligned}
```

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

function survival(m::Weibull,age) 
    return exp(-cumhazard(m,age))
end

"""
    InverseWeibull(;m,σ)

Construct a mortality model following Weibull's law of mortality.

The Inverse-Weibull proves useful for modelling the childhood and teenage years, because the logarithm of h(x) is a concave function.
 - `m >0` is a measure of location
 - `σ >0` is measure of dispersion

```math
\\begin{aligned}
\\mathrm{hazard} \\left( {\\rm age} \\right) &= \\frac{\\frac{1}{\\sigma} \\cdot \\left( \\frac{{\\rm age}}{m} \\right)^{\\frac{ - m}{\\sigma} - 1}}{e^{\\left( \\frac{{\\rm age}}{m} \\right)^{\\frac{ - m}{\\sigma}}} - 1}
\\\\
\\mathrm{cumhazard}\\left( {\\rm age} \\right) &=  - \\log\\left( 1 - e^{ - \\left( \\frac{{\\rm age}}{m} \\right)^{\\frac{ - m}{\\sigma}}} \\right)
\\\\
\\mathrm{survival}\\left( {\\rm age} \\right) &=  e^{ - \\mathrm{cumhazard}\\left( m, {\\rm age} \\right)}
\\end{aligned}
```

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

function survival(m::InverseWeibull,age) 
    return exp(-cumhazard(m,age))
end

"""
    Perks(;a,b,c,d)

Construct a mortality model following Perks' law of mortality.

``
\\mathrm{hazard} \\left( {\\rm age} \\right) = \\frac{a + b \\cdot c^{{\\rm age}}}{b \\cdot c^{ - {\\rm age}} + 1 + d \\cdot c^{{\\rm age}}}
``

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
    VanderMaen(;a,b,c,i,n)

Construct a mortality model following VanderMaen's law of mortality.

``
\\mathrm{hazard} \\left( {\\rm age} \\right) = a + b \\cdot {\\rm age} + c \\cdot {\\rm age}^{2} + \\frac{i}{n - {\\rm age}}
``

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
    VanderMaen2(;a,b,i,n)

Construct a mortality model following VanderMaen2's law of mortality.

``
\\mathrm{hazard} \\left( {\\rm age} \\right) = a + b \\cdot {\\rm age} + \\frac{i}{n - {\\rm age}}
``

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
    StrehlerMildvan(;k,v₀,b,d)

Construct a mortality model following StrehlerMildvan's law of mortality.

``
\\mathrm{hazard} \\left( {\\rm age} \\right) = k \\cdot e^{\\frac{\\left(  - v_0 \\right) \\cdot \\left( 1 - b \\cdot {\\rm age} \\right)}{d}}
``

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
    Beard(;a,b,k)

Construct a mortality model following Beard's law of mortality.

``
\\mathrm{hazard} \\left( {\\rm age} \\right) = \\frac{a \\cdot e^{b \\cdot {\\rm age}}}{1 + k \\cdot a \\cdot e^{b \\cdot {\\rm age}}}
``

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
    MakehamBeard(;a,b,c,k)

Construct a mortality model following MakehamBeard's law of mortality.

``
\\mathrm{hazard} \\left( {\\rm age} \\right) =\\left( {\\rm age} \\right) = \\frac{a \\cdot e^{b \\cdot {\\rm age}}}{1 + k \\cdot a \\cdot e^{b \\cdot {\\rm age}}} + c
``

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
    Quadratic(;a,b,c)

Construct a mortality model following Quadratic law of mortality.

``
\\mathrm{hazard} \\left( {\\rm age} \\right) = a + b \\cdot {\\rm age} + c \\cdot {\\rm age}^{2}
``

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
    GammaGompertz(;a,b,γ)

Construct a mortality model following GammaGompertz law of mortality.

``
\\mathrm{hazard} \\left( {\\rm age} \\right) = \\frac{a \\cdot e^{b \\cdot {\\rm age}}}{1 + \\frac{a \\cdot \\gamma}{b} \\cdot \\left( e^{b \\cdot {\\rm age}} - 1 \\right)}
``

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
    Siler(;a,b,c,d,e)

Construct a mortality model following Siler law of mortality.

``
\\mathrm{hazard} \\left( {\\rm age} \\right) = a \\cdot e^{\\left(  - b \\right) \\cdot {\\rm age}} + c + d \\cdot e^{e \\cdot {\\rm age}}
``

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
    HeligmanPollard(;a,b,c,d,e,f,g,h)

Construct a mortality model following HeligmanPollard law of mortality with 8 parameters.

``
\\mathrm{hazard} \\left( {\\rm age} \\right) = a \\cdot e^{\\left(  - b \\right) \\cdot {\\rm age}} + c + d \\cdot e^{e \\cdot {\\rm age}}
``


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
    HeligmanPollard2(;a,b,c,d,e,f,g,h)

Construct a mortality model following HeligmanPollard (alternate) law of mortality with 8 parameters.


```math
\\begin{aligned}
\\mu_1 &= a^{\\left( {\\rm age} + b \\right)^{c}} + \\frac{g \\cdot h^{{\\rm age}}}{1 + g \\cdot h^{{\\rm age}}}
\\\\
\\mu_2 &= d \\cdot e^{\\left(  - e \\right) \\cdot \\left( \\log\\left( \\frac{{\\rm age}}{f} \\right) \\right)^{2}}
\\\\
\\mathrm{hazard}\\left( {\\rm age} \\right) &= \\begin{cases}
\\mu_1 & \\text{if } \\left( {\\rm age} = 0 \\right)\\\\
\\mu_1 + \\mu_2 & \\text{otherwise}
\\end{cases}
\\end{aligned}
```

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
    HeligmanPollard3(;a,b,c,d,e,f,g,h,k)

Construct a mortality model following HeligmanPollard (alternate) law of mortality with 9 parameters.

```math
\\begin{aligned}
\\mu_1 &= a^{\\left( {\\rm age} + b \\right)^{c}} + \\frac{g \\cdot h^{{\\rm age}}}{1 + k \\cdot g \\cdot h^{{\\rm age}}}
\\\\
\\mu_2 &= d \\cdot e^{\\left(  - e \\right) \\cdot \\left( \\log\\left( \\frac{{\\rm age}}{f} \\right) \\right)^{2}}
\\\\
\\mathrm{hazard}\\left( {\\rm age} \\right) &= \\begin{cases}
\\mu_1 & \\text{if } \\left( {\\rm age} = 0 \\right)\\\\
\\mu_1 + \\mu_2 & \\text{otherwise}
\\end{cases}
\\end{aligned}
```

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
    HeligmanPollard4(;a,b,c,d,e,f,g,h,k)

Construct a mortality model following HeligmanPollard (alternate) law of mortality with 9 parameters.

```math
\\begin{aligned}
\\mu_1 &= a^{\\left( {\\rm age} + b \\right)^{c}} + \\frac{g \\cdot h^{{\\rm age}^{k}}}{1 + g \\cdot h^{{\\rm age}^{k}}}
\\\\
\\mu_2 &= d \\cdot e^{\\left(  - e \\right) \\cdot \\left( \\log\\left( \\frac{{\\rm age}}{f} \\right) \\right)^{2}}
\\\\
\\mathrm{hazard}\\left( {\\rm age} \\right) &= \\begin{cases}
\\mu_1 & \\text{if } \\left( {\\rm age} = 0 \\right)\\\\
\\mu_1 + \\mu_2 & \\text{otherwise}
\\end{cases}
\\end{aligned}
```

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
    RogersPlanck(;a₀, a₁, a₂, a₃, a, b, c, d, u)

Construct a mortality model following RogersPlanck law of mortality.

``
\\mathrm{hazard}\\left( {\\rm age} \\right) = a_0 + a_1 \\cdot e^{\\left(  - a \\right) \\cdot {\\rm age}} + a_2 \\cdot e^{b \\cdot \\left( {\\rm age} - u \\right) - e^{\\left(  - c \\right) \\cdot \\left( {\\rm age} - u \\right)}} + a_3 \\cdot e^{d \\cdot {\\rm age}}
``

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
    Martinelle(;a,b,c,d,k)

Construct a mortality model following Martinelle's law of mortality.

``
\\mathrm{hazard}\\left( {\\rm age} \\right) = \\frac{a \\cdot e^{b \\cdot {\\rm age}} + c}{1 + d \\cdot e^{b \\cdot {\\rm age}}} + k \\cdot e^{b \\cdot {\\rm age}}
``

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
    Kostaki(;a,b,c,d,e1,e2,f,g,h)

Construct a mortality model following Kostaki's law of mortality. A nine-parameter adaptation of `HeligmanPollard`.

```math
\\begin{aligned}
\\mu_1 &= a^{\\left( {\\rm age} + b \\right)^{c}} + g \\cdot h^{{\\rm age}}
\\\\
\\mu_2 &= \\begin{cases}
d \\cdot e^{ - \\left( e1 \\cdot \\log\\left( \\frac{{\\rm age}}{f} \\right) \\right)^{2}} & \\text{if } \\left( {\\rm age} \\leq f \\right)\\\\
d \\cdot e^{ - \\left( e2 \\cdot \\log\\left( \\frac{{\\rm age}}{f} \\right) \\right)^{2}} & \\text{otherwise}
\\end{cases}
\\\\
\\eta &= \\begin{cases}
\\mu_1 & \\text{if } \\left( {\\rm age} = 0 \\right)\\\\
\\mu_1 + \\mu_2 & \\text{otherwise}
\\end{cases}
\\\\
\\mathrm{hazard}\\left( {\\rm age} \\right) &= \\frac{\\eta}{1 + \\eta}

\\end{aligned}
```

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


> Kostaki, A. (1992). A nine‐parameter version of the Heligman‐Pollard formula. Mathematical Population Studies, 3(4), 277–288. doi:10.1080/08898489209525346 
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
    Kannisto(;a,b)

Construct a mortality model following Kannisto's law of mortality.

```math
\\begin{aligned}
\\mathrm{hazard}\\left( {\\rm age} \\right) &= \\frac{a \\cdot e^{b \\cdot {\\rm age}}}{1 + a \\cdot e^{b \\cdot {\\rm age}}}
\\\\
\\mathrm{cumhazard}\\left( {\\rm age} \\right) &= 1/a * log((1 + b*exp(b*age)) / (1 + a))
\\\\
\\mathrm{survival}\\left( {\\rm age} \\right) &= e^{ - \\mathrm{cumhazard}\\left( m, {\\rm age} \\right)}
\\end{aligned}
```

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

function  survival(m::Kannisto,age)
    return exp(-cumhazard(m,age))
end


"""
    KannistoMakeham(;a,b,c)

Construct a mortality model following KannistoMakeham's law of mortality.

``
\\mathrm{hazard}\\left( {\\rm age} \\right) = \\frac{a \\cdot e^{b \\cdot {\\rm age}}}{1 + a \\cdot e^{b \\cdot {\\rm age}}} + c
``

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

"""
    Carriere1(;p₁,p₂,μ₁,μ₂,μ₃,σ₁,σ₂,σ₃)
Construct a mortality model following Carriere1's law of mortality.
``
\\mathrm{hazard}\\left( {\\rm age} \\right) = 
``
Default args:
    p₁ = 0.003
    p₂ = 0.007
    μ₁ = 2.7
    μ₂ = 3
    μ₃ = 88
    σ₁ = 15
    σ₂ = 6
    σ₃ = 9.5
"""
Base.@kwdef struct Carriere1 <: ParametricMortality
    p₁ = 0.003
    p₂ = 0.007
    μ₁ = 2.7
    μ₂ = 3
    μ₃ = 88
    σ₁ = 15
    σ₂ = 6
    σ₃ = 9.5
end

function hazard(m::Carriere1,age)
    c = cumhazard(m, age)
    return vcat(c[1], diff(c))
end

function cumhazard(m::Carriere1,age)
    return -log.(survival(m, age))
end

function survival(m::Carriere1,age)
    @unpack p₁, p₂, μ₁, μ₂, μ₃, σ₁, σ₂, σ₃ = m

    s_weibull = exp.(-(age ./ μ₁) .^ (μ₁ / σ₁))
    s_inv_weibull = 1.0 .- exp.(-(age ./ μ₂) .^ (-μ₂ / σ₁))
    s_alter_gompertz = survival(AlterGompertz(μ=μ₃, σ=σ₃), age)

    f1 = max(0.0001, min(p₁, 1.0))
    f2 = max(0.0001, min(p₂, 1.0))
    f3 = 1.0 - f1 - f2

    return max.(0.0, min.(1.0, f1 .* s_weibull .+ f2 .* s_inv_weibull .+ f3 .* s_alter_gompertz))
end

"""
    Carriere2(;p₁,p₂,μ₁,μ₂,μ₃,σ₁,σ₂,σ₃)
Construct a mortality model following Carriere2's law of mortality.
``
\\mathrm{hazard}\\left( {\\rm age} \\right) = 
``
Default args:
    p₁ = 0.01
    p₂ = 0.01
    μ₁ = 1
    μ₂ = 49
    μ₃ = 49
    σ₁ = 2
    σ₂ = 7
    σ₃ = 7
"""
Base.@kwdef struct Carriere2 <: ParametricMortality
    p₁ = 0.01
    p₂ = 0.01
    μ₁ = 1
    μ₂ = 49
    μ₃ = 49
    σ₁ = 2
    σ₂ = 7
    σ₃ = 7
end

function hazard(m::Carriere2,age)
    c = cumhazard(m, age)
    return vcat(c[1], diff(c))
end

function cumhazard(m::Carriere2,age)
    return -log.(survival(m, age))
end

function survival(m::Carriere2,age)
    @unpack p₁, p₂, μ₁, μ₂, μ₃, σ₁, σ₂, σ₃ = m

    s_weibull = exp.(-(age ./ μ₁) .^ (μ₁ / σ₁))
    s_inv_gompertz = survival(InverseGompertz(μ=μ₂, σ=σ₂), age)
    s_alter_gompertz = survival(AlterGompertz(μ=μ₃, σ=σ₃), age)

    f1 = max(0.0001, min(p₁, 1.0))
    f2 = max(0.0001, min(p₂, 1.0))
    f3 = 1.0 - f1 - f2

    return max.(0.0, min.(1.0, f1 .* s_weibull .+ f2 .* s_inv_gompertz .+ f3 .* s_alter_gompertz))
end

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
decrement(m::ParametricMortality,from_age,to_age) = 1 - survival(m, from_age, to_age)
decrement(m::ParametricMortality,to_age) = 1 - survival(m, to_age)

(m::ParametricMortality)(x) = μ(m, x)
Base.getindex(m::ParametricMortality,x) = m(x)
Base.broadcastable(pm::ParametricMortality) = Ref(pm)
