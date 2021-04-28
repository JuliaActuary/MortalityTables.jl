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

### Bundled [mort.SOA.org](https://mort.soa.org) Tables

Comes with some tables built in via [mort.SOA.org](https://mort.soa.org) and by using [you agree to their terms](https://mort.soa.org/TermsOfUse.aspx).

!!! note
    Many mortality tables from `mort.SOA.org` have been tested to work, however not *all* mortality tables have been tested. Additionally, `mort.SOA.org` provides non-mortality rate tables which may not be propoerly parsed by this package. Try it and look at a few of the resulting values. If you encounter any issues, please [file an issue](https://github.com/JuliaActuary/MortalityTables.jl/issues).

Sample of some of the included table sets:

```plaintext
2017 Loaded CSO
2015 VBT
2001 VBT
2001 CSO
1980 CSO
1980 CET
```

[Click here to see the full list of tables included.](https://github.com/JuliaActuary/MortalityTables.jl/blob/master/BundledTables.md)

If you would like more tables added by default, please open a GitHub issue with the request.

### Other [mort.SOA.org](https://mort.soa.org) Tables

Given a table id ([for example](https://mort.soa.org/ViewTable.aspx?&TableIdentity=60029) `60029`)
you can request the table directly from the SOA's mortality table service. Remember
that not all tables have been tested, though the standard source format should mean
compatibility with `MortalityTables.jl`.

```julia
aus_life_table_female = get_SOA_table(60029)
aus_life_table_female[0]  # returns the attained age 0 rate of 0.10139
```

You can combine it with the bundled tables too:

```julia
tables = MortalityTables.tables()

get_SOA_table!(tables,60029) # this modifies `tables` by adding the new table

t = tables["Australian Life Tables 1891-1900 Female"]
t[0]  # returns the attained age 0 rate of 0.10139
```

### Load custom set of tables

Download the `.xml` [aka the (`XTbML` format)](https://mort.soa.org/About.aspx) version of the table from [mort.SOA.org](https://mort.soa.org) and place it in a directory of your choosing. Then call `MortaliyTables.tables(path_to_your_dir)`.


### From CSV

If you have a CSV file that is from [mort.SOA.org](https://mort.SOA.org), or follows the same structure, then you can load and parse the table into a `MortalityTable` like so:

```julia
using CSV
using MortalityTables

path = "path/to/table.csv"
file = CSV.File(path,header=false) # SOA's CSV files have no true header
table = MortalityTable(file)
```

### From XTbML

If you have a file using the XTbML format:

```
using MortalityTables
path = "path/to/table.xml"
table = MortalityTables.readXTbML(path)
```

### Custom Tables

Say you have an ultimate vector and select matrix, and you want to leverage the MortalityTables package.

Here's an example, where we first construct the `UlitmateMortality` and then combine
it with the select rates to get a `SelectMortality` table.

```julia
using MortalityTables

# represents attained ages of 15 through 100
ult_vec = [0.005, 0.008, ...,0.805,1.00]
ult = UltimateMortality(ult_vec,start_age = 15)
```

We can now use this the ulitmate rates all by itself:

```julia
q(ult,15,1) # 0.005
```

And join with the select rates, which for our example will start at age 0:

```julia
# attained age going down the column, duration across
select_matrix = [ 0.001 0.002 ... 0.010;
                  0.002 0.003 ... 0.012;
                  ...
                ]
sel_start_age = 0
sel = SelectMortality(select_matrix,ult,start_age = 0)

sel[0][0] #issue age 0, attained age 0 rate of  0.001
sel[0][100] #issue age 0, attained age 100 rate of  1.0
```

Lastly, to take the `SelectMortality` and `UltimateMortality` we just created,
we can combine them into one stored object, along with a `TableMetaData`:

```julia
my_table = MortalityTable(
              s1,
              u1,
              metadata=TableMetaData(name="My Table", comments="Rates for Product XYZ")
              )
```

```@docs
MortalityTables.readXTbML
MortalityTables.table
MortalityTables.get_SOA_table
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

## Rates, Survival and Decrement

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

### Docstrings

```@docs
MortalityTables.survival
MortalityTables.decrement
```

## Life Expectancy

Calculate curtate or complete life expectancy.

### Docstrings

```@docs
MortalityTables.life_expectancy
```

## Fractional Year Assumptions

When evaluating survival over partial years when you are given full year mortality
rates, you must make an assumption over how those deaths are distributed throughout
the year. Three assumptions are provided as options and are based on formulas
from the [2016 Experience Study Calculations paper from the SOA](https://www.soa.org/globalassets/assets/Files/Research/2016-10-experience-study-calculations.pdf), specifically pages 40-44.

The three assumptions are:

- `Uniform()` which assumes an increasing force of mortality throughout the year.
- `Constant()` which assumes a level force of mortality throughout the year.
- `Balducci()` which assumes a decreasing force of mortality over the year. It seems [to
be for making it easier](https://www.soa.org/globalassets/assets/library/research/actuarial-research-clearing-house/1978-89/1988/arch-1/arch88v17.pdf) to calculate successive months by hand rather than any theoretical basis.

```@docs
MortalityTables.DeathDistribution
MortalityTables.Balducci
MortalityTables.Uniform
MortalityTables.Constant
```