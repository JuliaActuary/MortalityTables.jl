# Parametric Mortality Models

```@meta
CurrentModule = MortalityTables
```

```@contents
Pages = ["ParametricMortalityModels.md"]
Depth = 5
```

## Introduction

Parametric mortality models are not collections of annual rates, but rather construct a mathematical representation of a mortality curve.

They are often used underlying traditional table construction, but are converted to a list of annual rates for convenience in practice. `MortalityTables.jl` provides them and are essentially interchangeable versus the traditional table structure in most use cases.

Many of these models were adapted from the [MortalityLaws](https://github.com/mpascariu/MortalityLaws) R package, by Marius Pascariu & Vladimir Canudas-Romo.

## Usage Example

Use like so:

```julia
a = 0.0002
b = 0.13
c = 0.001
m = MortalityTables.Makeham(a=a,b=b,c=c)
g = MortalityTables.Gompertz(a=a,b=b)
```

Now some examples with `m`, but could use `g` interchangeably:

```julia
age = 20
hazard(m,20)          # the force of mortality at age 20 (equivalently `m(20)`)
decrement(m,20,25)    # the five year cumulative mortality rate
survival(m,20,25) # the five year survival rate
```

### Other notes

- Because of the large number of models and the likelihood for overlap with other things (e.g. `Quadratic` or `Weibull` would be expected to be found in other contexts as well), these models are not exported from the package (except `Makeham` and `Gompertz`), so you need to call them by prefixing with `MortalityTables`.
  - e.g. : `MortalityTables.Kostaki()`
- Because of the large number of parameters for the models, the arguments are keyword rather than positional: `MortalityTables.Gompertz(a=0.01,b=0.2)`
- Keywords with Unicode names also accept ASCII spellings: `sigma` for `σ`, `gamma` for `γ`, `v0` for `v₀`, and `a0` to `a3` for `a₀` to `a₃`. Passing both spellings of one keyword is an `ArgumentError`.
- The alternative Makeham/Gompertz models also accept `mu` for `μ`; the Carrière models accept `p1`, `p2`, `mu1` to `mu3`, and `sigma1` to `sigma3`.
- The models have default values, so they can be called without args like this: `MortalityTables.Gompertz()`.
  - See the help text for what the default values are: `?Gompertz`

## Functions

Every model defines its `hazard`. `cumhazard`, `survival` and `decrement` follow from it, and MortalityTables integrates the hazard numerically for a model without a closed-form `cumhazard`.

```@docs; canonical=false
MortalityTables.hazard
MortalityTables.cumhazard
MortalityTables.μ
```

Also see [`survival`](@ref) and [`decrement`](@ref).

## Alternative parameterizations and survival mixtures

[`Makeham2`](@ref) and [`Gompertz2`](@ref) describe the Gompertz component using a location `μ` and a scale `σ`. They are equivalent to the existing laws with `a = exp(-μ/σ)/σ` and `b = 1/σ`. They retain the location and scale directly so that an underflowed coefficient `a` does not erase mortality near age `μ`.

```@example carriere
using MortalityTables
g = MortalityTables.Gompertz2(mu = 88, sigma = 9.5)
(hazard(g, 80), decrement(g, 80, 81))
```

[`Carriere1`](@ref) mixes Weibull, inverse-Weibull, and Gompertz survival functions. [`Carriere2`](@ref) substitutes inverse Gompertz for the middle component. The initial weights are `p₁`, `p₂`, and `1-p₁-p₂`; each component's survival at age zero is one. These are mixtures of populations with different lifetimes: the hazard at an attained age weights the component hazards by the population still surviving at that age.

```@example carriere
m = MortalityTables.Carriere1(p1 = 0.003, p2 = 0.007, sigma2 = 6)
ages = [20.0, 40.0, 60.0, 80.0]
(hazards = hazard.(m, ages), annual_decrements = decrement.(m, ages, ages .+ 1))
```

The force of mortality `hazard(m, age)` is instantaneous; `decrement(m, age, age + 1)` is the one-year probability of death conditional on survival to `age`. `cumhazard`, interval survival, and interval decrements use the same continuous law, including at fractional ages.

Carrière weights must be nonnegative and sum to at most one for the first two components. Zero weights are supported. Locations and scales must be finite and positive; invalid inputs raise an `ArgumentError`. A negative age raises a `DomainError`. There is no clipping of weights or survival probabilities.

At age zero, an active Weibull component with `μ₁ < σ₁` gives an infinite hazard, while survival is exactly one and cumulative hazard is zero. This singularity is integrable. The default `Carriere1` also has an infinite mean lifetime because its inverse-Weibull component has `μ₂/σ₂ ≤ 1`; use parameters with `μ₂/σ₂ > 1` or set `p₂ = 0` when a finite mean is needed. `Carriere2` has a finite mean for all allowed parameters.

The defaults are inherited from MortalityLaws and are examples rather than a calibrated mortality table. In particular, the default `Carriere1` hazard declines sharply at old ages as inverse-Weibull survivors dominate.

`Carriere1` uses `σ₂` for the inverse-Weibull component. Older MortalityLaws implementations used `σ₁` and ignored `σ₂`, before the [v2.3.0 development correction](https://github.com/mpascariu/MortalityLaws/commit/8f13abb1417cc22aec01b82848eea689bf01f048). Parameters fitted with that older formula generally need refitting. MortalityLaws' reported hazards are finite differences of cumulative hazard; `hazard` here is the instantaneous force.

The zero-age limits above apply to ordinary numeric ages. ForwardDiff age derivatives at exactly zero remain unsupported for `InverseWeibull` and mixtures containing it: they can return `NaN` even when the right-hand limit is finite.

## Available Models

These models are subtypes of `ParametricMortality`:

```@index
Pages   = ["ParametricMortalityModels.md"]
Modules = [MortalityTables]
Order   = [:type]
```

```@docs; canonical=false
MortalityTables.Makeham
MortalityTables.Gompertz
MortalityTables.Makeham2
MortalityTables.Gompertz2
MortalityTables.Carriere1
MortalityTables.Carriere2
MortalityTables.InverseGompertz
MortalityTables.Opperman
MortalityTables.Thiele
MortalityTables.Wittstein
MortalityTables.Weibull
MortalityTables.InverseWeibull
MortalityTables.Perks
MortalityTables.VanderMaen
MortalityTables.VanderMaen2
MortalityTables.StrehlerMildvan
MortalityTables.Beard
MortalityTables.MakehamBeard
MortalityTables.Quadratic
MortalityTables.GammaGompertz
MortalityTables.Siler
MortalityTables.HeligmanPollard
MortalityTables.HeligmanPollard2
MortalityTables.HeligmanPollard3
MortalityTables.HeligmanPollard4
MortalityTables.RogersPlanck
MortalityTables.Martinelle
MortalityTables.Kostaki
MortalityTables.Kannisto
MortalityTables.KannistoMakeham
```
