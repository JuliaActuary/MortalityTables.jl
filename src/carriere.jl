# Alternative Gompertz parameterization and Carrière survival mixtures.

"""
    Makeham2(; μ=49, σ=7.692308, c=0.001)

Makeham's law in location/scale parameters, with force of mortality
``h(x) = e^{(x-μ)/σ}/σ + c``. It is equivalent to
`Makeham(a=exp(-μ/σ)/σ, b=1/σ, c=c)`.

`μ` must be finite, `σ` finite and positive, and `c` finite and nonnegative.
`mu` and `sigma` are ASCII aliases for `μ` and `σ`; passing both spellings of a
keyword is an `ArgumentError`. Ages must be nonnegative.

The location and scale are retained as fields: transforming them to `a` and `b`
can underflow `a` even when mortality near age `μ` is representable.
"""
struct Makeham2{T<:Real} <: ParametricMortality
    μ::T
    σ::T
    c::T
    function Makeham2(μ::Real, σ::Real, c::Real)
        μ, σ, c = promote(μ, σ, c)
        isfinite(_primal(μ)) || throw(ArgumentError("μ must be finite"))
        isfinite(_primal(σ)) && _primal(σ) > 0 || throw(ArgumentError("σ must be finite and positive"))
        isfinite(_primal(c)) && _primal(c) >= 0 || throw(ArgumentError("c must be finite and nonnegative"))
        return new{typeof(μ)}(μ, σ, c)
    end
end
Makeham2(; μ = _Unset(), mu = _Unset(), σ = _Unset(), sigma = _Unset(), c = 0.001) =
    Makeham2(_keyword(:μ => μ, :mu => mu, 49), _keyword(:σ => σ, :sigma => sigma, 7.692308), c)

"""
    Gompertz2(; μ=49, σ=7.7)

Gompertz's law in location/scale parameters, with force of mortality
``h(x) = e^{(x-μ)/σ}/σ``. Returns a [`Makeham2`](@ref) with `c=0`, equivalent to
`Gompertz(a=exp(-μ/σ)/σ, b=1/σ)`. `mu` and `sigma` are ASCII keyword aliases.
"""
function Gompertz2(; μ = _Unset(), mu = _Unset(), σ = _Unset(), sigma = _Unset())
    location = _keyword(:μ => μ, :mu => mu, 49)
    scale = _keyword(:σ => σ, :sigma => sigma, 7.7)
    return Makeham2(location, scale, zero(scale))
end

_check_nonnegative_age(age) = _primal(age) >= 0 ? nothing : throw(DomainError(age, "age must be nonnegative"))

function hazard(m::Makeham2, age)
    _check_nonnegative_age(age)
    return exp((age - m.μ) / m.σ - log(m.σ)) + m.c
end

function cumhazard(m::Makeham2, from, to)
    _check_nonnegative_age(from)
    _check_nonnegative_age(to)
    from == to == Inf && return zero(float(from + to + m.c))
    to < from && return -cumhazard(m, to, from)
    Δ = to - from
    t = Δ / m.σ
    y = (to - m.μ) / m.σ
    rate = exp(y - log(m.σ))
    # Integrate the endpoint force over an effective interval width. This is
    # smooth at Δ = 0, even when exp(y) overflows but exp(y)/σ is finite.
    # For long intervals the width tends to σ; forming Δ/t could give Inf/Inf.
    growth = if isfinite(rate) &&
                (floatmin(typeof(_primal(rate))) <= _primal(rate) || iszero(_primal(Δ)))
        width = t < 1 ? Δ * _exprel(-t) : m.σ * -expm1(-t)
        rate * width
    elseif abs(t) < 1 && (from - m.μ) / m.σ < 1
        # Keep a subnormal force from losing digits before a large width rescales
        # it, while preserving the local derivatives in the small-growth tail.
        exp((from - m.μ) / m.σ) * t * _exprel(t)
    elseif iszero(Δ)
        zero(t)
    else
        # The force itself can under/overflow while its integral is representable.
        exp(y + log(-expm1(-t)))
    end
    return growth + _times_age(m.c, Δ)
end
cumhazard(m::Makeham2, age) = cumhazard(m, zero(age), age)

# Validate positional and keyword construction alike; never clip a probability or a scale.
function _check_carriere_parameters(p₁, p₂, μ₁, μ₂, μ₃, σ₁, σ₂, σ₃)
    p₁, p₂ = _primal(p₁), _primal(p₂)
    isfinite(p₁) && isfinite(p₂) && p₁ >= 0 && p₂ >= 0 && p₁ + p₂ <= 1 ||
        throw(ArgumentError("weights must satisfy p₁ ≥ 0, p₂ ≥ 0, and p₁ + p₂ ≤ 1"))
    all(x -> isfinite(_primal(x)) && _primal(x) > 0, (μ₁, μ₂, μ₃, σ₁, σ₂, σ₃)) ||
        throw(ArgumentError("locations and scales must be finite and positive"))
    return nothing
end

"""
    Carriere1(; p₁=0.003, p₂=0.007, μ₁=2.7, μ₂=3, μ₃=88, σ₁=15, σ₂=6, σ₃=9.5)

Carrière's three-component survival mixture: Weibull, inverse Weibull, and Gompertz.
For weights ``w = (p₁, p₂, 1-p₁-p₂)`` and component survival functions ``S_i``,

```math
S(x) = \\sum_{i=1}^3 w_i S_i(x), \\qquad
H(x) = -\\log S(x), \\qquad
h(x) = \\frac{\\sum_{i=1}^3 w_i S_i(x) h_i(x)}{S(x)}.
```

The components are `Weibull(m=μ₁, σ=σ₁)`, `InverseWeibull(m=μ₂, σ=σ₂)`, and
`Gompertz2(μ=μ₃, σ=σ₃)`. Each is normalized to survival one at age zero.
The hazard is an instantaneous force, not an annual mortality probability;
use `decrement(model, age, age + 1)` for the latter.

Weights must be nonnegative with `p₁ + p₂ ≤ 1`; zero weights are allowed.
All locations and scales must be finite and positive. Ages must be nonnegative.
The hazard at zero is infinite when `p₁ > 0` and `μ₁ < σ₁` (an integrable
Weibull singularity); survival there is one. The inverse-Weibull tail gives an
infinite mean lifetime when `p₂ > 0` and `μ₂/σ₂ ≤ 1`, as with the defaults.
The defaults are inherited from MortalityLaws, not a calibrated mortality table:
their hazard declines sharply at old ages as inverse-Weibull survivors dominate.

This implementation uses `σ₂` for the inverse-Weibull component, following
Carrière's paper. Older MortalityLaws implementations used `σ₁` there and ignored
`σ₂`; that wiring was corrected in MortalityLaws' v2.3.0 development commit
[`8f13abb`](https://github.com/mpascariu/MortalityLaws/commit/8f13abb1417cc22aec01b82848eea689bf01f048).
Parameters fitted with the older formula generally need refitting before use here.
MortalityLaws also reports finite differences of cumulative hazard, whereas
`hazard` here is the instantaneous force.

The zero-age limits above apply to ordinary numeric ages. ForwardDiff age
derivatives at exactly zero are not reliable for `InverseWeibull` or the Carrière
mixtures. They can return `NaN` where the right-hand derivative is finite, or
where it is `+Inf` because of a singular Weibull component (as in both mixtures'
defaults, including `Carriere2`, which has no inverse-Weibull component).

ASCII aliases are `p1`, `p2`, `mu1` to `mu3`, and `sigma1` to `sigma3`.
Passing both spellings of a keyword is an `ArgumentError`.

See Carrière (1992), [*Parametric Models for Life Tables*](https://www.soa.org/globalassets/assets/library/research/transactions-of-society-of-actuaries/1990-95/1992/january/tsa92v446.pdf),
Transactions of the Society of Actuaries 44, 77–99, and the
[MortalityLaws implementation](https://github.com/mpascariu/MortalityLaws).
"""
struct Carriere1{T<:Real} <: ParametricMortality
    p₁::T
    p₂::T
    μ₁::T
    μ₂::T
    μ₃::T
    σ₁::T
    σ₂::T
    σ₃::T
    function Carriere1(p₁::Real, p₂::Real, μ₁::Real, μ₂::Real, μ₃::Real, σ₁::Real, σ₂::Real, σ₃::Real)
        args = promote(p₁, p₂, μ₁, μ₂, μ₃, σ₁, σ₂, σ₃)
        _check_carriere_parameters(args...)
        return new{typeof(args[1])}(args...)
    end
end

"""
    Carriere2(; p₁=0.01, p₂=0.01, μ₁=1, μ₂=49, μ₃=49, σ₁=2, σ₂=7, σ₃=7)

Carrière's survival mixture of `Weibull(m=μ₁, σ=σ₁)`,
`InverseGompertz(m=μ₂, σ=σ₂)`, and `Gompertz2(μ=μ₃, σ=σ₃)`.
With ``w = (p₁,p₂,1-p₁-p₂)``, ``S = ∑ᵢ wᵢ Sᵢ``, ``H = -log(S)``, and
``h = (∑ᵢ wᵢ Sᵢ hᵢ)/S``. Each component is normalized at age zero.

Parameter validation, ASCII keyword aliases, and the hazard's limit at age zero
are as for [`Carriere1`](@ref). Unlike `Carriere1`, this mixture has a finite
mean lifetime for every allowed parameter set.
"""
struct Carriere2{T<:Real} <: ParametricMortality
    p₁::T
    p₂::T
    μ₁::T
    μ₂::T
    μ₃::T
    σ₁::T
    σ₂::T
    σ₃::T
    function Carriere2(p₁::Real, p₂::Real, μ₁::Real, μ₂::Real, μ₃::Real, σ₁::Real, σ₂::Real, σ₃::Real)
        args = promote(p₁, p₂, μ₁, μ₂, μ₃, σ₁, σ₂, σ₃)
        _check_carriere_parameters(args...)
        return new{typeof(args[1])}(args...)
    end
end

function _carriere_keywords(defaults;
    p₁ = _Unset(), p1 = _Unset(), p₂ = _Unset(), p2 = _Unset(),
    μ₁ = _Unset(), mu1 = _Unset(), μ₂ = _Unset(), mu2 = _Unset(), μ₃ = _Unset(), mu3 = _Unset(),
    σ₁ = _Unset(), sigma1 = _Unset(), σ₂ = _Unset(), sigma2 = _Unset(), σ₃ = _Unset(), sigma3 = _Unset())
    return (_keyword(:p₁ => p₁, :p1 => p1, defaults[1]),
            _keyword(:p₂ => p₂, :p2 => p2, defaults[2]),
            _keyword(:μ₁ => μ₁, :mu1 => mu1, defaults[3]),
            _keyword(:μ₂ => μ₂, :mu2 => mu2, defaults[4]),
            _keyword(:μ₃ => μ₃, :mu3 => mu3, defaults[5]),
            _keyword(:σ₁ => σ₁, :sigma1 => sigma1, defaults[6]),
            _keyword(:σ₂ => σ₂, :sigma2 => sigma2, defaults[7]),
            _keyword(:σ₃ => σ₃, :sigma3 => sigma3, defaults[8]))
end
Carriere1(; kwargs...) = Carriere1(_carriere_keywords((0.003, 0.007, 2.7, 3, 88, 15, 6, 9.5); kwargs...)...)
Carriere2(; kwargs...) = Carriere2(_carriere_keywords((0.01, 0.01, 1, 49, 49, 2, 7, 7); kwargs...)...)

_carriere_components(m::Carriere1) =
    (Weibull(m = m.μ₁, σ = m.σ₁), InverseWeibull(m = m.μ₂, σ = m.σ₂), Gompertz2(μ = m.μ₃, σ = m.σ₃))
_carriere_components(m::Carriere2) =
    (Weibull(m = m.μ₁, σ = m.σ₁), InverseGompertz(m = m.μ₂, σ = m.σ₂), Gompertz2(μ = m.μ₃, σ = m.σ₃))
_carriere_weights(m) = (m.p₁, m.p₂, 1 - (m.p₁ + m.p₂))

# Testing for a weight equal to one alone is insufficient: 1 - a tiny positive
# weight can round to one while that tiny component still contributes mortality.
_sole_component(weights) = count(!iszero, weights) == 1 ? findfirst(!iszero, weights) : nothing

function _mixture_cumhazards(components, weights, age)
    # Ignore inactive components, including singularities and exhausted survival.
    return map(components, weights) do component, weight
        iszero(weight) ? oftype(float(age + weight), Inf) : cumhazard(component, age)
    end
end

function _mixture_cumhazard(H, weights)
    h0 = minimum(H)
    isinf(h0) && return h0
    # Factor out the largest component survival. Near zero, sum decrements instead
    # of subtracting an almost-one survival; in the tail, sum the scaled survivals.
    loss = sum(map((w, h) -> w * -expm1(h0 - h), weights, H))
    remainder = loss < 1 / 2 ? -log1p(-loss) : -log(sum(map((w, h) -> w * exp(h0 - h), weights, H)))
    return h0 + remainder
end

function cumhazard(m::Union{Carriere1,Carriere2}, age)
    _check_nonnegative_age(age)
    weights = _carriere_weights(m)
    H = _mixture_cumhazards(_carriere_components(m), weights, age)
    return _mixture_cumhazard(H, weights)
end

# Only used when *every* active H has overflowed at a finite age. Comparing log H
# still identifies the surviving component; any distinguishable larger H has
# vanishing relative survival. The inverse laws are then in their asymptotic tails.
_overflow_logH(m::Weibull, x) = (m.m / m.σ) * log(x / m.m)
_overflow_logH(m::InverseWeibull, x) = log(m.m) - log(m.σ) + log(log(x / m.m))
_overflow_logH(m::InverseGompertz, x) = log(x - m.m) - log(m.σ)
_overflow_logH(m::Makeham2, x) = (x - m.μ) / m.σ + log(-expm1(-x / m.σ))

function _mixture_masses(components, weights, age)
    H = _mixture_cumhazards(components, weights, age)
    h0 = minimum(H)
    if isinf(h0)
        logH = map((m, w) -> iszero(w) ? h0 : _overflow_logH(m, age), components, weights)
        smallest = minimum(logH)
        return map((w, h) -> h == smallest ? w : zero(w), weights, logH)
    end
    return map((w, h) -> w * exp(h0 - h), weights, H)
end

# Interval hazards after the starting cumulative hazard has overflowed. For Weibull,
# factor the difference of powers before exponentiating; the inverse tails have
# already reached their limiting power/exponential forms at floating-point precision.
function _overflow_interval(m::Weibull, from, to)
    k = m.m / m.σ
    t = k * log1p((to - from) / from)
    log_difference = t < 1 ? log(expm1(t)) : t + log1p(-exp(-t))
    return exp(k * log(from / m.m) + log_difference)
end
_overflow_interval(m::InverseWeibull, from, to) = (m.m / m.σ) * log1p((to - from) / from)
_overflow_interval(m::InverseGompertz, from, to) = (to - from) / m.σ
_overflow_interval(m::Makeham2, from, to) = cumhazard(m, from, to)

function cumhazard(m::Union{Carriere1,Carriere2}, from, to)
    _check_nonnegative_age(from)
    _check_nonnegative_age(to)
    (iszero(to - from) || from == to == Inf) && return zero(float(from + to + m.p₁))
    to < from && return -cumhazard(m, to, from)
    weights = _carriere_weights(m)
    components = _carriere_components(m)
    if isinf(cumhazard(m, from))
        masses = _mixture_masses(components, weights, from)
        posterior = map(w -> w / sum(masses), masses)
        H = map(components, posterior) do component, weight
            iszero(weight) ? oftype(float(from + weight), Inf) : _overflow_interval(component, from, to)
        end
        return _mixture_cumhazard(H, posterior)
    end
    i = _sole_component(weights)
    # Retain the component's interval kernel at the vertices of the simplex,
    # including Gompertz intervals where both endpoint cumulative hazards overflow.
    i === nothing || return cumhazard(components[i], from, to)
    return cumhazard(m, to) - cumhazard(m, from)
end

function hazard(m::Union{Carriere1,Carriere2}, age)
    _check_nonnegative_age(age)
    weights = _carriere_weights(m)
    components = _carriere_components(m)
    i = _sole_component(weights)
    i === nothing || return hazard(components[i], age)
    if isinf(age)
        # Among these components the slowest-decaying survival also has the
        # smallest limiting hazard (0, an exponential rate, or Inf).
        return minimum(map(components, weights) do component, weight
            iszero(weight) ? oftype(float(age + weight), Inf) : hazard(component, age)
        end)
    end
    masses = _mixture_masses(components, weights, age)
    total = sum(masses)
    return sum(map(components, masses) do component, mass
        iszero(mass) ? zero(total) : (mass / total) * hazard(component, age)
    end)
end
