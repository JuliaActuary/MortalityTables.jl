"""
    curtate_life_expectancy(rates, age)
    curtate_life_expectancy(law, age; rtol, atol)

The curtate expectation of life at `age`: the expected number of whole years lived after `age`,
``e_x = \\sum_{k \\ge 1} {}_kp_x``, where ``{}_kp_x`` is `survival(m, age, age + k)`.

For a vector of rates, `age` is a whole age in the table. The sum runs to the table's last age,
where the expectancy is zero, even if the table does not end with a rate of one. An age outside
the table is a `BoundsError`.

For a parametric law, `age` is any real age up to `omega(law)`. An age past `omega` is a
`DomainError`.
- With a finite `omega`, the sum runs to ``k = \\lfloor \\omega - x \\rfloor``. Survival left at
  `omega` counts as death there, so the expectancy at `omega` is zero.
- With `omega = Inf`, the first `K` terms are summed and the rest is estimated from the integral
  of survival ``S(t) = {}_tp_x``. Survival does not increase (the hazard is not negative), so the
  rest lies between ``\\int_{K+1}^\\infty S`` and ``\\int_K^\\infty S``, an interval no wider than
  ``S(K)``. The estimate uses its midpoint, integrated with QuadGK. `K` grows until ``S(K)/2``
  plus QuadGK's error estimate is at most `max(atol, rtol * abs(estimate))`. This is an
  estimated accuracy, not a certified bound. If 100,000 terms are not enough, it throws an
  `ArgumentError`.

`rtol` defaults to `sqrt(eps(T))`, where `T` is the type of the survival values, and `atol` to
`0`. They apply only to a law with `omega = Inf`, and must be finite and nonnegative (either may be
zero), or an `ArgumentError` names the one that is not.

See also [`complete_life_expectancy`](@ref).

# Examples
```julia-repl
julia> qs = UltimateMortality([0.1, 0.3, 0.6, 1]);

julia> curtate_life_expectancy(qs, 0) # 0.9 + 0.9 * 0.7 + 0.9 * 0.7 * 0.4
1.782

julia> curtate_life_expectancy(qs, 3)
0.0
```
"""
function curtate_life_expectancy(table::AbstractArray, age::Integer)
    # sum of the survival probabilities to each later whole age, accumulated in a single pass
    # rather than recomputed from `age` each time. `table[age]` is read first, so an age outside
    # the table is a BoundsError rather than an empty sum. The sum starts from the rates' numeric
    # type, not from `p`, whose rate can be `missing` at the last age.
    p = one(_rate_type(table)) * (1 - table[age])
    s = zero(_rate_type(table))
    for a in (age + 1):lastindex(table)
        s += p
        p *= 1 - table[a]
    end
    return s
end

# the type of a law's survival values from `age`
_expectancy_type(m::ParametricMortality, age) = typeof(survival(m, age, age))

# An infinite, NaN or negative tolerance would accept any error estimate, or none.
function _check_expectancy_tolerances(rtol, atol)
    for (name, value) in ((:rtol, rtol), (:atol, atol))
        value isa Real && isfinite(value) && value >= 0 ||
            throw(ArgumentError("$name must be finite and nonnegative; got $(repr(value))"))
    end
    return nothing
end

# the most terms the sum for a law with omega = Inf adds before it gives up
const _CURTATE_MAX_TERMS = 100_000

function curtate_life_expectancy(
        m::ParametricMortality, age::Real;
        rtol = sqrt(eps(_expectancy_type(m, age))), atol = 0,
    )
    _check_expectancy_tolerances(rtol, atol)
    # past a finite omega the sum below would be empty and give zero, so the age is checked
    _check_age(m, age)
    T = _expectancy_type(m, age)
    S(t) = survival(m, age, age + t)
    ω = omega(m)
    isfinite(ω) && return sum(S, 1:floor(Int, ω - age); init = zero(T))
    # s = Σ_{k=1}^{K} S(k) and SK = S(K), with K doubled until the midpoint estimate of the rest
    # is accurate enough
    s, SK, K = zero(T), one(T), 0
    while K < _CURTATE_MAX_TERMS
        K_next = min(max(2K, 1), _CURTATE_MAX_TERMS)
        for k in (K + 1):K_next
            SK = S(k)
            s += SK
        end
        K = K_next
        # the rest lies between ∫_{K+1}^∞ S and ∫_K^∞ S; their midpoint in one integral
        Q, E = quadgk(t -> (S(t) + S(t + 1)) / 2, T(K), T(Inf); rtol = rtol / 2, atol = atol / 2)
        estimate = s + Q
        isfinite(estimate) && isfinite(E) && SK / 2 + E <= max(atol, rtol * abs(estimate)) && return estimate
    end
    throw(ArgumentError("the curtate life expectancy did not converge within $(_CURTATE_MAX_TERMS) terms"))
end

"""
    complete_life_expectancy(rates, age, dist = UniformDeaths())
    complete_life_expectancy(law, age, dist = UniformDeaths(); rtol, atol, maxevals)

The complete expectation of life at `age`: the expected remaining lifetime,
``\\mathring{e}_x = \\int_0^{\\omega - x} {}_tp_x \\, dt``, where ``{}_tp_x`` is
`survival(m, age, age + t)`.

For a vector of rates, survival within each year of age follows the `DeathDistribution` `dist`,
and survival is integrated one year of age at a time. A fractional `age` conditions on survival
to `age`, as `survival` does. The integral runs to the table's last age, where the expectancy is
zero. An age outside the table is a `BoundsError`.

For a parametric law, survival is integrated with QuadGK up to `omega(law)`. `dist` is accepted
and ignored, since a law is continuous. With a finite `omega`, survival left at `omega` counts as
death there, so the expectancy at `omega` is zero. An age past `omega` is a `DomainError`.
`rtol`, `atol` and `maxevals` go to QuadGK, with its defaults: `atol = 0`, `rtol` the square root
of the machine epsilon of the integration interval's type when `atol` is zero, and
`maxevals = 10^7`. The tolerances must be finite and nonnegative and `maxevals` a positive integer,
or an `ArgumentError` names the one that is not. The integral is returned only if it and QuadGK's
error estimate are finite and the estimate is at most `max(atol, rtol * abs(integral))`; the
estimate is QuadGK's, not a proof of accuracy. Otherwise an `ArgumentError` reports that the
integration did not converge: the expectancy may be infinite, as for a law whose survival falls too
slowly, or the integral may need a larger `maxevals`.

See also [`curtate_life_expectancy`](@ref).

# Examples
```julia-repl
julia> qs = UltimateMortality([0.1, 0.3, 0.6, 1]);

julia> complete_life_expectancy(qs, 0) # 0.95 + 0.9 * 0.85 + 0.9 * 0.7 * 0.7
2.156

julia> complete_life_expectancy(qs, 0.5, ConstantForce())
1.7198453593325438
```
"""
function complete_life_expectancy(table::AbstractArray, age::Real, dist::DeathDistribution = UniformDeaths())
    # past the last age the loop below is empty, and the year containing an age just past it
    # has a rate, so the age is checked against the table
    firstindex(table) <= age <= lastindex(table) || throw(BoundsError(table, age))
    # the integral of survival from `age` over each year of age, in which survival is smooth,
    # times `p`, the survival from `age` to the start of that year
    T = _survival_type(table, age, float(age))
    e, p = zero(T), one(T)
    for to in (floor(Int, age) + 1):lastindex(table)
        from = max(age, oftype(age, to - 1))
        e += p * quadgk(u -> survival(table, from, u, dist), from, to)[1]
        p *= survival(table, from, to, dist)
    end
    return e
end

function complete_life_expectancy(
        m::ParametricMortality, age::Real, ::DeathDistribution = UniformDeaths();
        atol = 0, rtol = iszero(atol) ? sqrt(eps(float(typeof(omega(m) - age)))) : 0, maxevals = 10^7,
    )
    _check_expectancy_tolerances(rtol, atol)
    maxevals isa Integer && maxevals > 0 || throw(ArgumentError("maxevals must be a positive integer; got $(repr(maxevals))"))
    # Integrate over the remaining lifetime t. Integrating over ages [age, omega] gives the same bits
    # but is about 10% slower. The default tolerances are QuadGK's, for the interval's type.
    I, E = quadgk(t -> survival(m, age, age + t), 0, omega(m) - age; atol, rtol, maxevals)
    tol = max(atol, rtol * abs(I))
    isfinite(I) && isfinite(E) && E <= tol && return I
    throw(ArgumentError("the complete life expectancy did not converge: integral $I, estimated error $E, tolerance $tol"))
end
