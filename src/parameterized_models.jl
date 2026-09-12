# Reference/Adapted from:
# https://mortality.org/File/GetDocument/Public/HMD_4th_Symposium/Pascariu_poster.pdf
# https://github.com/mpascariu/MortalityLaws/blob/master/R/MortalityLaw_models.R

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
struct Makeham{T<:Real} <: ParametricMortality
    a::T
    b::T
    c::T
end
Makeham(; a=0.0002, b=0.13, c=0.001) = Makeham(promote(a, b, c)...)

function hazard(m::Makeham,age)
    (; a, b, c) = m
    return a*exp(b*age) + c
end

function cumhazard(m::Makeham,age)
    (; a, b, c) = m
    return a / b * (exp(b*age) - 1) + age * c
end


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
struct InverseGompertz{T<:Real} <: ParametricMortality
    m::T
    σ::T
end
InverseGompertz(; m=49, σ=7.7) = InverseGompertz(promote(m, σ)...)


function hazard(model::InverseGompertz,age)
    (; m, σ) = model
    return 1 / σ * exp(-(age - m)/σ) / (exp(exp(-(age - m)/σ)) - 1)
end

function cumhazard(model::InverseGompertz,age)
    (; m, σ) = model
    # negative log of the closed-form survival function
    return -log((1 - exp(-exp(-(age - m)/σ))) / (1 - exp(-exp(m/σ))))
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
struct Opperman{T<:Real} <: ParametricMortality
    a::T
    b::T
    c::T
end
Opperman(; a=0.04, b=0.0004, c=0.001) = Opperman(promote(a, b, c)...)

function hazard(m::Opperman,age) 
    (; a, b, c) = m
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
struct Thiele{T<:Real} <: ParametricMortality
    a::T
    b::T
    c::T
    d::T
    e::T
    f::T
    g::T
end
Thiele(; a=0.02474, b=0.3, c=0.004, d=0.5, e=25, f=0.0001, g=0.13) = Thiele(promote(a, b, c, d, e, f, g)...)

function hazard(m::Thiele,age) 
    (; a, b, c, d, e, f, g) = m
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
struct Wittstein{T<:Real} <: ParametricMortality
    a::T
    b::T
    n::T
    m::T
end
Wittstein(; a=1.5, b=1., n=0.5, m=100) = Wittstein(promote(a, b, n, m)...)

function hazard(model::Wittstein,age)
    (; a, b, m, n) = model
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
struct Weibull{T<:Real} <: ParametricMortality
    m::T
    σ::T
end
Weibull(; m=1.0, σ=2.0) = Weibull(promote(m, σ)...)

function hazard(model::Weibull,age)
    (; m, σ) = model
    if age == 0
        return 1.0
    else 
        return 1 / σ * (age / m)^(m / σ - 1)
    end
end

function cumhazard(model::Weibull,age)
    (; m, σ) = model
    return (age / m) ^ (m / σ)
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
struct InverseWeibull{T<:Real} <: ParametricMortality
    m::T
    σ::T
end
InverseWeibull(; m=5.0, σ=10.0) = InverseWeibull(promote(m, σ)...)

function hazard(model::InverseWeibull,age)
    (; m, σ) = model
    return (1/σ) * (age/m)^(-m/σ - 1) / (exp((age/m)^(-m/σ)) - 1)
end

function cumhazard(model::InverseWeibull,age)
    (; m, σ) = model
    return -log(1 - exp(-(age/m)^(-m/σ)))
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
struct Perks{T<:Real} <: ParametricMortality
    a::T
    b::T
    c::T
    d::T
end
Perks(; a=0.002, b=0.13, c=0.01, d=0.01) = Perks(promote(a, b, c, d)...)

function hazard(m::Perks,age) 
    (; a, b, c, d) = m
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
struct VanderMaen{T<:Real} <: ParametricMortality
    a::T
    b::T
    c::T
    i::T
    n::T
end
VanderMaen(; a=0.01, b=1., c=0.01, i=100., n=200.) = VanderMaen(promote(a, b, c, i, n)...)

function hazard(m::VanderMaen,age)
    (; a, b, c, i, n) = m
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
struct VanderMaen2{T<:Real} <: ParametricMortality
    a::T
    b::T
    i::T
    n::T
end
VanderMaen2(; a=0.01, b=1., i=100., n=200.) = VanderMaen2(promote(a, b, i, n)...)

function hazard(m::VanderMaen2,age)
    (; a, b, i, n) = m
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
struct StrehlerMildvan{T<:Real} <: ParametricMortality
    k::T
    v₀::T
    b::T
    d::T
end
StrehlerMildvan(; k=0.01, v₀=2.5, b=0.2, d=6.0) = StrehlerMildvan(promote(k, v₀, b, d)...)

function hazard(m::StrehlerMildvan,age)
    (; k, v₀, b, d) = m
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
struct Beard{T<:Real} <: ParametricMortality
    a::T
    b::T
    k::T
end
Beard(; a=0.002, b=0.13, k=1.) = Beard(promote(a, b, k)...)

function hazard(m::Beard,age)
    (; a, b, k) = m
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
struct MakehamBeard{T<:Real} <: ParametricMortality
    a::T
    b::T
    c::T
    k::T
end
MakehamBeard(; a=0.002, b=0.13, c=0.01, k=1.) = MakehamBeard(promote(a, b, c, k)...)

function hazard(m::MakehamBeard,age)
    (; a, b, c, k) = m
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
struct Quadratic{T<:Real} <: ParametricMortality
    a::T
    b::T
    c::T
end
Quadratic(; a=0.01, b=1., c=0.01) = Quadratic(promote(a, b, c)...)

function hazard(m::Quadratic,age)
    (; a, b, c) = m
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
struct GammaGompertz{T<:Real} <: ParametricMortality
    a::T
    b::T
    γ::T
end
GammaGompertz(; a=0.002, b=0.13, γ=1) = GammaGompertz(promote(a, b, γ)...)

function hazard(m::GammaGompertz,age)
    (; a, b, γ) = m
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
struct Siler{T<:Real} <: ParametricMortality
    a::T
    b::T
    c::T
    d::T
    e::T
end
Siler(; a=0.0002, b=0.13, c=0.001, d=0.001, e=0.013) = Siler(promote(a, b, c, d, e)...)

function hazard(m::Siler,age)
    (; a, b, c, d, e) = m
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
struct HeligmanPollard{T<:Real} <: ParametricMortality
    a::T
    b::T
    c::T
    d::T
    e::T
    f::T
    g::T
    h::T
end
HeligmanPollard(; a=.0005, b=.004, c=.08, d=.001, e=10, f=17, g=.00005, h=1.1) = HeligmanPollard(promote(a, b, c, d, e, f, g, h)...)

function hazard(m::HeligmanPollard,age)
    (; a, b, c, d, e, f, g, h) = m
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
struct HeligmanPollard2{T<:Real} <: ParametricMortality
    a::T
    b::T
    c::T
    d::T
    e::T
    f::T
    g::T
    h::T
end
HeligmanPollard2(; a=0.0005, b=0.004, c=0.08, d=0.001, e=10., f=17., g=0.00005, h=1.1) = HeligmanPollard2(promote(a, b, c, d, e, f, g, h)...)

function hazard(m::HeligmanPollard2,age)
    (; a, b, c, d, e, f, g, h) = m
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
struct HeligmanPollard3{T<:Real} <: ParametricMortality
    a::T
    b::T
    c::T
    d::T
    e::T
    f::T
    g::T
    h::T
    k::T
end
HeligmanPollard3(; a=.0005, b=.004, c=.08, d=.001, e=10, f=17, g=.00005, h=1.1, k=1.) = HeligmanPollard3(promote(a, b, c, d, e, f, g, h, k)...)

function hazard(m::HeligmanPollard3,age)
    (; a, b, c, d, e, f, g, h, k) = m
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
struct HeligmanPollard4{T<:Real} <: ParametricMortality
    a::T
    b::T
    c::T
    d::T
    e::T
    f::T
    g::T
    h::T
    k::T
end
HeligmanPollard4(; a=.0005, b=.004, c=.08, d=.001, e=10, f=17, g=.00005, h=1.1, k=1.) = HeligmanPollard4(promote(a, b, c, d, e, f, g, h, k)...)

function hazard(m::HeligmanPollard4,age)
    (; a, b, c, d, e, f, g, h, k) = m
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
struct RogersPlanck{T<:Real} <: ParametricMortality
    a₀::T
    a₁::T
    a₂::T
    a₃::T
    a::T
    b::T
    c::T
    d::T
    u::T
end
RogersPlanck(; a₀=0.0001, a₁=0.02, a₂=0.001, a₃=0.0001, a=2., b=0.001, c=100., d=0.1, u=0.33) = RogersPlanck(promote(a₀, a₁, a₂, a₃, a, b, c, d, u)...)

function hazard(m::RogersPlanck,age) 
    (; a₀, a₁, a₂, a₃, a, b, c, d, u) = m
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
struct Martinelle{T<:Real} <: ParametricMortality
    a::T
    b::T
    c::T
    d::T
    k::T
end
Martinelle(; a=0.001, b=0.13, c=0.001, d=0.1, k=0.001) = Martinelle(promote(a, b, c, d, k)...)

function hazard(m::Martinelle,age)
    (; a, b, c, d, k) = m
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
struct Kostaki{T<:Real} <: ParametricMortality
    a::T
    b::T
    c::T
    d::T
    e1::T
    e2::T
    f::T
    g::T
    h::T
end
Kostaki(; a=0.0005, b=0.01, c=0.10, d=0.001, e1=3., e2=0.1, f=25., g=.00005, h=1.1) = Kostaki(promote(a, b, c, d, e1, e2, f, g, h)...)

function hazard(m::Kostaki,age) 
    (; a, b, c, d, e1, e2, f, g, h) = m
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
\\mathrm{cumhazard}\\left( {\\rm age} \\right) &= \\frac{1}{b} \\log\\left( \\frac{1 + a \\cdot e^{b \\cdot {\\rm age}}}{1 + a} \\right)
\\\\
\\mathrm{survival}\\left( {\\rm age} \\right) &= e^{ - \\mathrm{cumhazard}\\left( m, {\\rm age} \\right)}
\\end{aligned}
```

Default args:

    a = 0.5
    b = 0.13
"""
struct Kannisto{T<:Real} <: ParametricMortality
    a::T
    b::T
end
Kannisto(; a=0.5, b=0.13) = Kannisto(promote(a, b)...)

function hazard(m::Kannisto,age)
    (; a, b) = m
    return  a * exp(b * age) / (1 + a * exp(b*age))
end

function cumhazard(m::Kannisto,age)
    (; a, b) = m
    return  1/b * log((1 + a*exp(b*age)) / (1 + a))
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
struct KannistoMakeham{T<:Real} <: ParametricMortality
    a::T
    b::T
    c::T
end
KannistoMakeham(; a=0.5, b=0.13, c=0.001) = KannistoMakeham(promote(a, b, c)...)

function hazard(m::KannistoMakeham,age)
    (; a, b, c) = m
    return  a * exp(b * age) / (1 + a * exp(b*age)) + c
end
