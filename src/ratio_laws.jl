# Reference/Adapted from:
# https://mortality.org/File/GetDocument/Public/HMD_4th_Symposium/Pascariu_poster.pdf
# https://github.com/mpascariu/MortalityLaws/blob/master/R/MortalityLaw_models.R

# The laws whose hazard is a bounded ratio of exponentials: Beard, Makeham–Beard, Gamma–Gompertz,
# Martinelle, Kannisto and Kannisto–Makeham. They share the ratio kernels in `parametric_kernels.jl`.

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
    return _logistic_ratio(a, k, _times_age(b, age))
end

"""
    MakehamBeard(;a,b,c,k)

Construct a mortality model following MakehamBeard's law of mortality.

``
\\mathrm{hazard} \\left( {\\rm age} \\right) = \\frac{a \\cdot e^{b \\cdot {\\rm age}}}{1 + k \\cdot a \\cdot e^{b \\cdot {\\rm age}}} + c
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
    return _logistic_ratio(a, k, _times_age(b, age)) + c
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
    iszero(a) && return zero(a * age)   # no hazard, even at an infinite age
    x = _times_age(b, age)
    # a·exp(x) / (1 + a·γ/b·expm1(x)), evaluated as a ratio in a (see `_amplitude_ratio`). Near
    # b·age = 0 it is a / (exp(-x) + a·γ·age·expm1(-x)/(-x)), with `_exprel`, which also gives the
    # b = 0 limit a / (1 + a·γ·age). Away from it the algebra follows the sign of x: for x > 0 it
    # is divided through by exp(x), a / (exp(-x) + a·γ/b·(1 - exp(-x))), finite where exp(x)
    # overflows (tending to b/γ); for x < 0 it is exp(x)·a / (1 + a·γ/b·expm1(x)), since there
    # exp(-x) can overflow (and a·γ/b·expm1(-x) would be 0·Inf when γ = 0).
    abs(x) < 1 && return _amplitude_ratio(a, _times_age(γ, age) * _exprel(-x), exp(-x))
    x > 0 && return _amplitude_ratio(a, -γ / b * expm1(-x), exp(-x))
    return exp(x) * _amplitude_ratio(a, γ / b * expm1(x), one(a))
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
    x = _times_age(b, age)
    # (a·exp(x) + c) / (1 + d·exp(x)) is bounded (it tends to a/d); for x > 0 it is divided
    # through by exp(x) so that it stays finite where exp(x) overflows
    bounded = x > 0 ? (a + c * exp(-x)) / (exp(-x) + d) : (a * exp(x) + c) / (1 + d * exp(x))
    return bounded + (iszero(k) ? zero(k * x) : k * exp(x))
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
    return _logistic_ratio(a, one(a), _times_age(b, age))
end

function cumhazard(m::Kannisto,age)
    (; a, b) = m
    iszero(a) && return zero(a * age)   # no hazard, even at an infinite age
    x = _times_age(b, age)
    # log((1 + a·exp(b·age)) / (1 + a)) / b. Away from b·age = 0, log(1 + v) with v = a·exp(x)
    # is log1p(v) where v ≤ 1, which is smooth in a at a = 0, and log1pexp(log(a) + x) where
    # v > 1, which cannot overflow and keeps the derivatives free of the large v. Near b·age = 0
    # it is log1p(u) / b with u = a·expm1(b·age) / (1 + a), which is a·age/(1 + a) at b = 0.
    if abs(x) >= 1
        v = a * exp(x)
        return ((v <= 1 ? log1p(v) : _log1pexp(log(a) + x)) - log1p(a)) / b
    end
    u = a / (1 + a) * expm1(x)
    return a / (1 + a) * age * _exprel(x) * _log1pdivx(u)
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
    return _logistic_ratio(a, one(a), _times_age(b, age)) + c
end
