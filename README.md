# MortalityTables

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaActuary.github.io/MortalityTables.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaActuary.github.io/MortalityTables.jl/dev)
![CI](https://github.com/JuliaActuary/MortalityTables.jl/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/JuliaActuary/MortalityTables.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaActuary/MortalityTables.jl)

A Julia package for working with MortalityTables. Has:
- First-class support for missing values.
- Lots of bundled SOA mort.soa.org tables
- Many common functions, including partial year mortality calculations (Uniform, Constant, Balducci)
- Friendly syntax and flexible usage

## Examples
### Quickstart

```julia
using MortalityTables

tables = MortalityTables.tables() # loads the tables stored in the package folder
vbt2001 = tables["2001 VBT Residual Standard Select and Ultimate - Male Nonsmoker, ANB"]

# indexed by issue age and duration for select rate
q(vbt2001.select,35,1)        # .00036
q(vbt2001.ultimate,95,1)          # .24298

# when accessing ultimate rates, don't always need to specify the duration
q(vbt2001.ultimate,95)          # .24298
```

### Example: Quickly access and compare tables
```julia
using MortalityTables, Plots


tables = MortalityTables.tables()
cso_2001 = tables["2001 CSO Super Preferred Select and Ultimate - Male Nonsmoker, ANB"]
cso_2017 = tables["2017 Loaded CSO Preferred Structure Nonsmoker Super Preferred Male ANB"]

issue_age = 27
durations = 1:30
mort = [q(cso_2001.select,issue_age,durations),
        q(cso_1980.ultimate,issue_age,durations)]
plot(
   durations,
   mort,
   label = ["2001 CSO M SuperPref NS" "1980 CSO M NS"],
   title = "Comparison of 1980 and 2001 CSO \n for 27-year-old male",
   xlabel="duration")
```
![Comparison of 2001 and 2017 CSO](https://user-images.githubusercontent.com/711879/80447339-c48e5180-88de-11ea-81ee-4babd0f84755.png)

Easily extend the analysis to move up the [ladder of abstraction](http://worrydream.com/LadderOfAbstraction/):

```julia
ages = 25:80
durs = 1:35

# compute the relative rates with the element-wise division ("brodcasting" in Julia)
diff = q(cso_2017.select,ages,durs) ./ q(cso_2001.select,ages,durs)

contour(durs,ages,diff,
        xlabel="duration",ylabel="issue ages",
        title="Relative difference between 2017 and 2001 CSO \n M PFN",
        fill=true
        )
```

![Heatmap comparison](https://user-images.githubusercontent.com/711879/80447494-251d8e80-88df-11ea-9335-761d1a1739c7.png)


## Usage

### Indexing

Tables are indexed by a starting age and duration (even ultimate tables, under the hood). For tables with a starting age that is defined, but you've requested an attained age beyond the defined rates, you will get a `BoundsError`. If you ask for a starting age that is not defined, you will get a `missing`.

The rationale for this is, for example, this [2001 CSO table](https://mort.soa.org/ViewTable.aspx?&TableIdentity=1076) is not defined for ages 15 and under and ends at 120.
- You will get a `missing` if you ask for starting age 10 or 15, because it's plausible that you could encounter a a starting age not defined by a table.
- You will get a `BoundsError` if you ask for an attained age 150 for someone select at age 16, because it is beyond the table's definition of its end.

#### Index by issue age and duration

```julia
q(vbt2001.select,35,1)        # .00036
q(vbt2001.select,35,61)       # .24298

# can easily get ranges of values:
q(vbt2001.select,35,1:30)     # [0.0036, 0.0048, ...]
```

#### Index by just age to get the ultimate rates
```julia
q(vbt2001,95)          # .24298
q(vbt2001,50:70)       # [0.00319, 0.00345, ...]
```

### Other Usage

#### Table Attributes
```julia
issue_age = 50
ω(vbt2001.select,issue_age)              # 120
omega(vbt2001.ultimate, issue_age)          # 120
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

#### Fractional Years
When evaluating survival over partial years when you are given full year mortality
rates, you must make an assumption over how those deaths are distributed throughout
the year. Three assumptions are provided as options and are based on formulas
from the [2016 Experience Study Calculations paper from the SOA](https://www.soa.org/globalassets/assets/Files/Research/2016-10-experience-study-calculations.pdf), specifically pages 40-44.

The three assumptions are:

- `Uniform()` which assumes an increasing force of mortality throughout the year.
- `Constant()` which assumes a level force of mortality throughout the year.
- `Balducci()` which assumes a decreasing force of mortality over the year. It seems [to
be for](https://www.soa.org/globalassets/assets/library/research/actuarial-research-clearing-house/1978-89/1988/arch-1/arch88v17.pdf) making it easier to calculate successive months by hand.

##### Usage

When you call a method below that uses the `time` argument (ie a period over which
    you want to calculate the force of mortality), if you use a non `Int` (integer)
    number then you need to specify the assumption as the last argument.

For example:

**_Don't_** need to specify because you gave an `Int` time:

```julia
# calculate the 5-year survival for a person issued at age 50 and in
# the first duration
issue_age = 50
duration = 1
time = 5
p(table,issue_age,duration,time)
```
**_Do_** need to specify because you gave an fractional floating time:

```julia
# calculate the 5-and-a-half-year survival for a
# person issued at age 50 and in the first duration
issue_age = 50
duration = 1
time = 5.5
p(table,issue_age,duration,time, Balducci())
```

Note that if you are passing floating point numbers as the `time`
argument, you still have to specify a mortality assumption because
[floating point numbers](https://en.wikipedia.org/wiki/Floating-point_arithmetic) are often *technically* not whole numbers even
if you define a number to be, say, `5.0`.


### Actuarial Notation Equivalants
#### `ₜp₍ₓ₎₊ₛ`
The probability that a life aged `x + s` who was select
at age `x` survives to at least age `x+s+t`.

```julia
issue_age = x
duration = s - 1
time = t
p(table,issue_age,duration,time)
```

#### `ₜpₓ`
The probability that a life aged `x` survives to at least age `t`.

```julia
issue_age = x
duration = 1
time = t
p(table,issue_age,duration,time)
```

#### `pₓ`
The probability that a life aged `x` survives through age `x+1`

```julia
issue_age = x
duration = 1
time = 1
p(table,issue_age,duration,time)
```

#### `ₜqₓ`
The probability that a life aged `x` dies by age `x+t`

```julia
issue_age = x
duration = 1
time = t
q(table,issue_age,duration,time)
```

#### `qₓ`
The probability that a life aged `x` dies by age `x+1`

```julia
issue_age = x
duration = 1
time = 1
q(table,issue_age,duration,time)
```

#### `ω`
Returns the last attained age which the table has defined (ie not including) for a given issue_age.

```julia
omega(table, issue_age)
ω(table, issue_age)
```

### Some Batteries Included

Comes with some tables built in via [mort.SOA.org](https://mort.soa.org) and by using [you agree to their terms](https://mort.soa.org/TermsOfUse.aspx).

Not all tables have been tested that they work by default, though I have not encountered issues with any of the the VBT/CSO/other usual tables.

Sample of some of the included table sets:
```
2017 Loaded CSO
2015 VBT
2001 VBT
2001 CSO
1980 CSO
1980 CET
```

[Click here to see the full list of tables included.](BundledTables.md)


## Adding more tables

### Getting tables from [mort.SOA.org](https://mort.soa.org)

Given a table id ([for example](https://mort.soa.org/ViewTable.aspx?&TableIdentity=60029) `60029`)
you can request the table directly from the SOA's mortality table service.

```
aus_life_table_female = get_SOA_table(60029)
q(aus_life_table_female,0,1)  # returns the issue age 0, first duration rate of 0.10139
```

You can combine it with the bundled tables too:

```
tables = MortalityTables.tables()

get_SOA_table!(tables,60029) # this modifies `tables` by adding the new table

t = tables["Australian Life Tables 1891-1900 Female"]
q(t,0,1)  # returns the issue age 0, first duration rate of 0.10139
```


### Constructing Dynamically

Say you have an ultimate vector and select matrix, and you want to leverage the MortalityTables package.

Here's an example, where we first construct the `UlitmateMortality` and then combine
it with the select rates to get a `SelectMortality` table.
```julia
using MortalityTables

# represents attained ages of 15 through 100
ult_start_age = 15
ult_vec = [0.005, 0.008, ...,0.805,1.00]
ult = UltimateMortality(ult_vec,ult_start_age)
```

We can now use this the ulitmate rates all by itself:
```julia
q(ult,15,1) # 0.005
q(ult,10,1) # missing (no age 10 rate)
```
And join with the select rates, which for our example will start at age 0:
```julia
# attained age going down the column, duration across
select_matrix = [ 0.001 0.002 ... 0.010;
                  0.002 0.003 ... 0.012;
                  ...
                ]
sel_start_age = 0
sel = SelectMortality(select_matrix,ult,sel_start_age)

q(sel,0,1) # 0.001
```

Lastly, to take the `SelectMortality` and `UltimateMortality` we just created,
we can combine them into one stored object, along with MetaData:

```julia
my_table = MortalityTable(
              s1,
              u1,
              TableMetaData(name="My Table", comments="Rates for Product XYZ")
              )
```



### Load with bundled tables

To add more tables for your use when loading with all of the other bundled tables, download the `.xml` [aka the (`XTbML` format)](https://mort.soa.org/About.aspx) version of the table from [mort.SOA.org](https://mort.soa.org) and place it in the directory the package is installed in. This is usually `~user/.julia/packages/MortalityTables/[changing hash value]/src/tables/`.

> :warning: *updating the package may remove your existing tables. Make a backup before updating your packages*

After placing packages in the folder above, restart Julia and the should be discoverable when you run `mt.Tables()`


### Getting more tables bundled with the package

If you would like more tables added by default, please open a GitHub issue with the request.


## References
- [SOA Mortality Tables](https://mort.soa.org/)
- [Actuarial Mathematics for Life Contingent Risks, 2nd ed](https://www.cambridge.org/vi/academic/subjects/statistics-probability/statistics-econometrics-finance-and-insurance/actuarial-mathematics-life-contingent-risks-2nd-edition?format=HB)
- [Experience Study Calculations, SOA](https://www.soa.org/globalassets/assets/Files/Research/2016-10-experience-study-calculations.pdf)

## Similar Projects
 - [Pyliferisk, a Python package](https://github.com/franciscogarate/pyliferisk)
