# Tables



```@meta
CurrentModule = MortalityTables
```

```@contents
Pages = ["Tables.md"]
Depth = 5
```

## Loading Tables

There are a variety of ways to load your own tables, reference online tables, or load bundled tables:

!!! Many mortality tables from `mort.SOA.org` have been tested to work, however not *all* mortality tables have been tested. Additionally, `mort.SOA.org` provides non-mortality rate tables which may not be propoerly parsed by this package. Try it and look at a few of the resulting values. If you encounter any issues, please [file an issue](https://github.com/JuliaActuary/MortalityTables.jl/issues).

```@docs
MortalityTables.readXTbML
MortalityTables.get_SOA_table
MortalityTables.tables
MortalityTables.get_SOA_table!
```

## Table Constructors

Use these to build your own `MortalityTables.jl`-compatible table:

```@docs
MortalityTables.MortalityTable
MortalityTables.SelectMortality
MortalityTables.UltimateMortality
MortalityTables.mortality_vector
```

## Table Attributes

Basic metadata about the table (automatically populated for some tables).

```@docs
MortalityTables.TableMetaData
```

Find the final age for which a table defines a rate.

```@docs
MortalityTables.omega
```

## Survival and Decrement and Rates

To access the rates, simply index by the attained age. Example:

```julia
julia> vbt2001 = tables["2001 VBT Residual Standard Select and Ultimate - Male Nonsmoker, ANB"]
MortalityTable (Insured Lives Mortality):
   Name:
       2001 VBT Residual Standard Select and Ultimate - Male Nonsmoker, ANB
   Fields:
       (:select, :ultimate, :metadata)
   Provider:
       Society of Actuaries
   mort.SOA.org ID:
       1118
   mort.SOA.org link:
       https://mort.soa.org/ViewTable.aspx?&TableIdentity=1118
   Description:
       2001 Valuation Basic Table (VBT) Residual Standard Select and Ultimate Table -  Male Nonsmoker.
       Basis: Age Nearest Birthday. 
       Minimum Select Age: 0. 
       Maximum Select Age: 99. 
       Minimum Ultimate Age: 25. 
       Maximum Ultimate Age: 120
```

The package revolves around easy-to-access vectors which are indexed by attained age:

```julia
julia> vbt2001.select[35]          # vector of rates for issue age 35
 0.00036
 0.00048
 ⋮
 0.94729
 1.0
 
julia> vbt2001.select[35][35]      # issue age 35, attained age 35
 0.00036
 
julia> vbt2001.select[35][50:end] # issue age 35, attained age 50 through end of table
0.00316
0.00345
 ⋮
0.94729
1.0

julia> vbt2001.ultimate[95]        # ultimate vectors only need to be called with the attained age
 0.24298
```

To calculate the survivorship or cumulative decrement:

```@docs
MortalityTables.survival
MortalityTables.decrement
```

## Fractional Year Assumptions

```@docs
MortalityTables.DeathDistribution
MortalityTables.Balducci
MortalityTables.Uniform
MortalityTables.Constant
```