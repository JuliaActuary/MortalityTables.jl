# MortalityTables

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://alecloudenback.github.io/MortalityTables.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://alecloudenback.github.io/MortalityTables.jl/dev)
[![Build Status](https://travis-ci.com/JuliaActuary/MortalityTables.jl.svg?branch=master)](https://travis-ci.com/JuliaActuary/MortalityTables.jl)
[![codecov](https://codecov.io/gh/JuliaActuary/MortalityTables.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaActuary/MortalityTables.jl)

A Julia package for working with MortalityTables. Has first-class support for missing values.

## Examples
### Quickstart

```julia
using MortalityTables

tables = MortalityTables.tables() # loads the tables stored in the package folder
vbt2001 = tables["2001 VBT Residual Standard Select and Ultimate - Male Nonsmoker, ANB"]

# indexed by issue age and duration for select rate
qx(vbt2001,35,1)        # .00036

# indexed by just attained age for ultimate rate
qx(vbt2001,95)          # .24298
```

### Example: Quickly access and compare tables
```julia
using MortalityTables, Plots


tables = MortalityTables.tables()
cso_2001 = tables["2001 CSO Super Preferred Select and Ultimate - Male Nonsmoker, ANB"]
cso_1980 = tables["1980 CSO - Male Nonsmoker, ANB"]

age = 27
durations = 1:30
plot(
   durations,
   [qx(cso_2001,age,durations),qx(cso_1980,age,durations)],
   label = ["2001 CSO M SuperPref NS" "1980 CSO M NS"],
   title = "Comparison of 1980 and 2001 CSO \n for 27-year-old male",
   xlabel="duration")
```
![Comparison of 2001 and 1980 CSO](https://user-images.githubusercontent.com/711879/74115814-789cfd80-4b76-11ea-9d27-c4901844bcb4.png)


## Usage

### Indexing

Tables are indexed by a starting age and duration (even ultimate tables, under the hood). For tables with a starting age that is defined, but you've requested an attained age beyond the defined rates, you will get a `BoundsError`. If you ask for a starting age that is not defined, you will get a `missing`.

The rationale for this is, for example, this [2001 CSO table](https://mort.soa.org/ViewTable.aspx?&TableIdentity=1076) is not defined for ages 15 and under.
- You will get a `missing` if you ask for starting age 10 or 150, because it's plausible that you could encounter a a starting age not defined by a table.
- You will get a `BoundsError` if you ask for an attained age 150 for someone select at age 16, because it is beyond the table's definition of its end.

#### Index by issue age and duration to get select rates

```julia
qx(vbt2001,35,1)        # .00036
qx(vbt2001,35,61)       # .24298

# can easily get ranges of values:
qx(vbt2001,35,1:30)     # [0.0036, 0.0048, ...]
```

#### Index by just age to get the ultimate rates
```julia
qx(vbt2001,95)          # .24298
qx(vbt2001,50:70)       # [0.00319, 0.00345, ...]
```

### Other Usage

#### Table Attributes
```julia
ω(vbt2001)              # 120
omega(vbt2001)          # 120
```

#### Table MetaData

When you have an expression that shows a Mortality table, it displays relevant information:

```julia
tables = MortalityTables.tables()
vbt2001 = tables["2001 VBT Residual Standard Select and Ultimate - Male Nonsmoker, ANB"]
```

This shows the following in a notebook or REPL:

```
MortalityTable:
   Name:
       2001 VBT Residual Standard Select and Ultimate - Male Nonsmoker, ANB
   Provider:
       Society of Actuaries
   mort.SOA.org ID:
       1118
   mort.SOA.org link:
       https://mort.soa.org/ViewTable.aspx?&TableIdentity=1118
   Description:
       2001 Valuation Basic Table (VBT) Residual Standard Select and Ultimate Table -  Male Nonsmoker. Basis: Age Nearest Birthday. Minimum Select Age: 0. Maximum Select Age: 99. Minimum Ultimate Age: 25. Maximum Ultimate Age: 120
```

### Exported functions
```julia
"""
ₜp₍ₓ₎₊ₛ , or the probability that a life aged `x + s` who was select
at age `x` survives to at least age `x+s+t`
"""
p(table::MortalityTable,x,s,t)


"""
ₜpₓ , or the probability that a life aged `x` survives to at least age `t`
"""
p(table::MortalityTable,x,t)


"""
pₓ , or the probability that a life aged `x` survives through age `x+1`
"""
p(table::MortalityTable,x)


"""
ₜq₍ₓ₎₊ₛ , or the probability that a life aged `x + s` who was select
at age `x` dies by least age `x+s+t`
"""
q(table::MortalityTable,x,s,t)


"""
ₜqₓ , or the probability that a life aged `x` dies by age `x+t`
"""
q(table::MortalityTable,x,t)


"""
qₓ , or the probability that a life aged `x` dies by age `x+1`
"""
q(table::MortalityTable,x)


"""
`qx` is a convenience function that allows you to get the rate at a given `age`.
If wanting select/ultimate rates, specify the `duration` and `age` should be the issue age.
"""
qx(table::MortalityTable,age)
qx(table::MortalityTable,age,duration)

"""
`omega` (also `ω`) returns the last attained age which the table has defined (ie not including)
`missing`
"""
omega(table::MortalityTable)
ω(table::MortalityTable)
```

### Some Batteries Included

Comes with some tables built in via [mort.SOA.org](https://mort.soa.org) and by using [you agree to their terms](https://mort.soa.org/TermsOfUse.aspx).

Not all tables have been tested that they work by default, though I have not encountered issues with any of the the VBT/CSO/other usual tables.

Included:
```
2017 Loaded CSO
2015 VBT
2001 VBT
2001 CSO
1980 CSO
1980 CET
```

[Click here to see the full list of tables included.](BundledTables.md)



### Adding more tables

To add more tables for your use, download the `.xml` (aka the (`Xtbml` format)[https://mort.soa.org/About.aspx]) version of the table from [mort.SOA.org](https://mort.soa.org) and place it in the directory the package is installed in. This is usually `~user/.julia/packages/MortalityTables/[changing hash value]/src/tables/`. *Note: updating the package may remove your existing tables. Make a backup before updating your packages*

After placing packages in the folder above, restart Julia and the should be discoverable when you run `mt.Tables()`

### Todos

- Docs
- Automatically parse built-in tables
- Add more built-in tables
- Usage Examples
- More tests
- Performance testing
- Easy assumption overlay (e.g. scalars)


### References
- [SOA Mortality Tables](https://mort.soa.org/)
- [Actuarial Mathematics for Life Contingent Risks, 2nd ed](https://www.cambridge.org/vi/academic/subjects/statistics-probability/statistics-econometrics-finance-and-insurance/actuarial-mathematics-life-contingent-risks-2nd-edition?format=HB)
- [Experience Study Calculations, SOA](https://www.soa.org/globalassets/assets/Files/Research/2016-10-experience-study-calculations.pdf)

### Similar Projects
 - [Pyliferisk, a Python package](https://github.com/franciscogarate/pyliferisk)
