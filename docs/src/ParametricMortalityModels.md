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
m[20]                 # the mortality rate at age 20
decrement(m,20,25)    # the five year cumulative mortality rate
survival(m,20,25) # the five year survival rate
```

### Other notes

- Because of the large number of models and the likelihood for overlap with other things (e.g. `Quadratic` or `Weibull` would be expected to be found in other contexts as well), these models Are not exported from the package, so you need to call them by prefixing with `MortalityTables`. 
  - e.g. : `MortalityTables.Kostaki()`
- Because of the large number of parameters for the models, the arguments are keyword rather than positional: `MortalityTables.Gompertz(a=0.01,b=0.2)`
- The models have default values, so they can be called without args like this: `MortalityTables.Gompertz()`.
  - See the help text for what the default values are: `?Gompertz`

## Functions

Many of the models only have analytical forms for either the `hazard` or `survival` curve. You can still call either on any of the models, as MortalityTables with numerically integrate or automatic differentiate as appropriate.

```@docs
decrement
MortalityTables.hazard
MortalityTables.cumhazard
MortalityTables.survival
MortalityTables.μ
```

## Available Models

These models are subtypes of `ParametricMortality`:

```@index
Pages   = ["ParametricMortalityModels.md"]
Modules = [MortalityTables]
Order   = [:type]
```

```@docs
MortalityTables.Makeham
MortalityTables.Gompertz
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