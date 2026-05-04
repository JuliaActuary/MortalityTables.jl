# Integration with FinanceCore.AbstractDeflator (FinanceCore 2.6+).
#
# Mortality processes are multiplicative-factor processes whose `factor` is the
# survival probability. Subtyping `FinanceCore.AbstractDeflator` lets every
# mortality type compose with yield curves, lapse forces, and default
# decrements via `FinanceCore.compose(...)`.
#
# Three concrete types live under `AbstractMortality`:
#
# 1. `MortalityVector` — wraps an age-indexed vector of `q`s. The canonical
#    return type from `UltimateMortality(...)` and `SelectMortality(...)[age]`.
#    Implements the `AbstractArray` interface by delegation so most existing
#    code that treated those constructors' returns as `OffsetArray` continues
#    to work.
#
# 2. `ParametricMortality` (Makeham, Gompertz, etc., defined in
#    parameterized_models.jl). Continuous-force closed-form models.
#
# 3. `ShiftedMortality` — wraps any `AbstractMortality` and shifts its age
#    origin to a valuation age. Constructed via `at_age(model, age)`. The
#    natural way to put a from-birth mortality model on a years-from-valuation
#    axis for composition with a yield curve.

"""
    AbstractMortality

Supertype for mortality models that act as multiplicative-factor processes
on the time axis.

`AbstractMortality <: FinanceCore.AbstractDeflator`, so every concrete
mortality model composes with yield curves, lapse forces, and default
decrements via `FinanceCore.compose(...)`.

Subtypes:
- [`MortalityVector`](@ref) — age-indexed vector of `q`s
- [`ParametricMortality`](@ref) — closed-form continuous models
- [`ShiftedMortality`](@ref) — origin-shifted wrapper for axis alignment
"""
abstract type AbstractMortality <: FinanceCore.AbstractDeflator end

# ─── MortalityVector ─────────────────────────────────────────────────────────

"""
    MortalityVector(qs::AbstractVector)

A mortality vector: an age-indexed vector of one-year decrement rates `qs`.
`MortalityVector` is the canonical return type from
[`UltimateMortality`](@ref) and from indexing a [`SelectMortality`](@ref) by
issue age.

`MortalityVector` implements the `AbstractArray` interface by delegating to
the wrapped vector, so `mv[age]`, `length(mv)`, iteration, and slicing all
work as if you had the underlying `OffsetArray` directly. `Base.parent(mv)`
returns the wrapped array for advanced use.

# Composition with a yield curve

`MortalityVector` indexing is on the **attained-age axis**: `mv[65]` is `q_65`,
`factor(mv, 65, 70)` is 5-year survival from age 65. To compose with a yield
curve on a **years-from-valuation axis**, either:

1. Use [`at_age`](@ref) to wrap the vector with a valuation-age origin:
   ```julia
   qs = UltimateMortality([…], start_age = 0)
   mort_at_65 = at_age(qs, 65)
   deflator = compose(yield_curve, mort_at_65)
   ```
2. Or use the explicit two-argument form `factor(mv, from_age, to_age)`
   throughout — no axis conversion needed.

The single-argument form `factor(mv, t)` anchors at the vector's `firstindex`,
matching the years-from-valuation convention when the vector starts at index 0.
"""
# MortalityVector subtypes AbstractMortality (and thus AbstractDeflator).
# We deliberately do NOT subtype AbstractArray so that the type remains
# rooted in the deflator hierarchy for `compose` semantics. AbstractArray-like
# operations (indexing, length, iteration, slicing) are forwarded by
# delegation methods, and we provide explicit `survival` methods that
# delegate to the existing AbstractArray implementations on the underlying
# storage.
struct MortalityVector{T, V <: AbstractVector{T}} <: AbstractMortality
    qs::V
end

Base.size(mv::MortalityVector) = size(mv.qs)
Base.length(mv::MortalityVector) = length(mv.qs)
Base.firstindex(mv::MortalityVector) = firstindex(mv.qs)
Base.lastindex(mv::MortalityVector) = lastindex(mv.qs)
Base.axes(mv::MortalityVector) = axes(mv.qs)
Base.eachindex(mv::MortalityVector) = eachindex(mv.qs)
Base.getindex(mv::MortalityVector, args...) = getindex(mv.qs, args...)
Base.iterate(mv::MortalityVector) = iterate(mv.qs)
Base.iterate(mv::MortalityVector, state) = iterate(mv.qs, state)
Base.parent(mv::MortalityVector) = mv.qs
Base.IteratorSize(::Type{<:MortalityVector}) = Base.HasLength()
Base.eltype(::Type{<:MortalityVector{T}}) where {T} = T

Base.:(==)(a::MortalityVector, b::MortalityVector) = a.qs == b.qs
Base.isequal(a::MortalityVector, b::MortalityVector) = isequal(a.qs, b.qs)

# Pretty-print: identify the type and the size, not the full vector
function Base.show(io::IO, mv::MortalityVector)
    print(io, "MortalityVector(<", length(mv.qs),
          " rates, ages ", firstindex(mv.qs), ":", lastindex(mv.qs), ">)")
end

# ─── factor, intensity, survival, decrement ──────────────────────────────────

# Two-arg factor: delegates to existing `survival(v, from, to)` on the
# underlying array.
function FinanceCore.factor(mv::MortalityVector, from_age::Integer, to_age::Integer)
    return survival(mv.qs, Int(from_age), Int(to_age))
end

# Two-arg factor with DeathDistribution: fractional ages.
function FinanceCore.factor(mv::MortalityVector, from_age::Real, to_age::Real, dd::DeathDistribution)
    return survival(mv.qs, from_age, to_age, dd)
end

# Two-arg factor with fractional ages and no DD raises informatively.
function FinanceCore.factor(mv::MortalityVector, from_age::Real, to_age::Real)
    throw(ArgumentError(
        "MortalityVector requires integer ages, or a DeathDistribution " *
        "(Uniform, Constant, Balducci) as the fourth argument for fractional ages."
    ))
end

# Single-arg factor: anchors at the vector's firstindex. For a vector with
# firstindex = 0 (the canonical years-from-valuation form, after `at_age`
# alignment or explicit construction), this gives "survival over t years
# from valuation age." Required for `pv(deflator, ::Cashflow)` to work.
function FinanceCore.factor(mv::MortalityVector, t::Real)
    isinteger(t) || throw(ArgumentError(
        "MortalityVector requires integer time arguments; got $t. " *
        "Pass a DeathDistribution as the third argument for fractional ages."))
    fi = firstindex(mv.qs)
    return survival(mv.qs, fi, fi + Int(t))
end

# survival on a MortalityVector is an alias for factor — keeps the actuarial
# idiom working without forcing users to swap to `factor`. We define explicit
# methods (rather than a vararg) to avoid ambiguity with the generic
# `survival(v, to_age, dd::DeathDistribution)` from MortalityTable.jl.
survival(mv::MortalityVector, from_age::Integer, to_age::Integer) =
    FinanceCore.factor(mv, from_age, to_age)
survival(mv::MortalityVector, from_age::Real, to_age::Real, dd::DeathDistribution) =
    FinanceCore.factor(mv, from_age, to_age, dd)
survival(mv::MortalityVector, from_age::Integer, to_age::Integer, dd::DeathDistribution) =
    FinanceCore.factor(mv, float(from_age), float(to_age), dd)

decrement(mv::MortalityVector, from_age::Real, to_age::Real) =
    1 - survival(mv, from_age, to_age)
decrement(mv::MortalityVector, from_age::Real, to_age::Real, dd::DeathDistribution) =
    1 - survival(mv, from_age, to_age, dd)

# omega: matches the existing convention (last attained age in the vector).
omega(mv::MortalityVector) = lastindex(mv.qs)

# ─── ShiftedMortality ────────────────────────────────────────────────────────

"""
    ShiftedMortality(model::AbstractMortality, age::Real)
    at_age(model::AbstractMortality, age::Real)

Wrap an [`AbstractMortality`](@ref) so that its origin is shifted to a
specified valuation `age`. The wrapper exposes a years-from-valuation axis:

```julia-repl
julia> mort = Makeham(a = 2.5e-5, b = 0.10, c = 1e-4);

julia> mort_at_65 = at_age(mort, 65);

julia> factor(mort_at_65, 5)        # 5-year survival from age 65
0.8973...
```

`at_age` is the idiomatic way to align a from-birth mortality model
([`ParametricMortality`](@ref)) or an absolute-age table ([`MortalityVector`](@ref))
to a years-from-valuation axis for composition with a yield curve:

```julia
deflator = compose(yield_curve, at_age(mort, 65))
```

`factor`, `intensity`, and `survival` on a `ShiftedMortality` delegate to
the wrapped model with shifted arguments.
"""
struct ShiftedMortality{M <: AbstractMortality, A <: Real} <: AbstractMortality
    model::M
    age::A
end

at_age(m::AbstractMortality, age::Real) = ShiftedMortality(m, age)

FinanceCore.factor(sm::ShiftedMortality, t) =
    FinanceCore.factor(sm.model, sm.age, sm.age + t)
FinanceCore.factor(sm::ShiftedMortality, from, to) =
    FinanceCore.factor(sm.model, sm.age + from, sm.age + to)
# intensity propagates only when the underlying model defines it.
FinanceCore.intensity(sm::ShiftedMortality, t) =
    FinanceCore.intensity(sm.model, sm.age + t)

survival(sm::ShiftedMortality, args...) = FinanceCore.factor(sm, args...)
decrement(sm::ShiftedMortality, args...) = 1 - FinanceCore.factor(sm, args...)

function Base.show(io::IO, sm::ShiftedMortality)
    print(io, "at_age(", repr(sm.model), ", ", sm.age, ")")
end

