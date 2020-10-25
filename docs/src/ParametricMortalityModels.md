# Parametric Mortality Models

Parametric mortality models are not collections of annual rates, but rather construct a mathematical representation of a mortality curve.

They are often used underlying traditional table construction, but are converted to a list of annual rates for convenience in practice. `MortalityTables.jl` provides them and are essentially interchangeable versus the traditional table structure in most use cases.

Many of these models were adapted from the [MortalityLaws](https://github.com/mpascariu/MortalityLaws) R package, by Marius Pascariu & Vladimir Canudas-Romo.

```@meta
CurrentModule = MortalityTables
```

```@contents
Pages = ["ParametricMortalityModels.md"]
Depth = 5
```

## Functions

```@docs
decrement
MortalityTables.hazard
MortalityTables.cumhazard
survival
MortalityTables.μ
```

## Available Models

These models are subtypes of `ParametricMortality`:

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