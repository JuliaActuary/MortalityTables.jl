# MortalityTables

#### Code Review: [![Build Status](https://travis-ci.org/alecloudenback/MortalityTables.jl.svg?branch=master)](https://travis-ci.org/alecloudenback/MortalityTables.jl) [![Coverage Status](https://coveralls.io/repos/github/alecloudenback/MortalityTables.jl/badge.svg?branch=master)](https://coveralls.io/github/alecloudenback/MortalityTables.jl?branch=master) [![codecov.io](http://codecov.io/github/alecloudenback/MortalityTables.jl/coverage.svg?branch=master)](http://codecov.io/github/alecloudenback/MortalityTables.jl?branch=master)
[![FOSSA Status](https://app.fossa.io/api/projects/git%2Bgithub.com%2Falecloudenback%2FMortalityTables.jl.svg?type=shield)](https://app.fossa.io/projects/git%2Bgithub.com%2Falecloudenback%2FMortalityTables.jl?ref=badge_shield)
A Julia package for working with MortalityTables. Has first-class support for missing values.

### Simple Usage Example

```julia
using MortalityTables
const mt = MortalityTables

tables = mt.Tables()
cso = tables["2001 CSO Super Preferred Select and Ultimate - Male Nonsmoker, ANB"]


mt.qx(cso,35,1) # == .00037
mt.qx(cso,35,61) # == .26719
mt.qx(cso,95) # == .26719
mt.qx(cso,35,95) # == missing (table doesn't extend that far)
```

### Five lines of code to visualize and compare two tables
*Well, eight if you count the import lines*
```julia
using Plots
using MortalityTables
const mt = MortalityTables

tbls = mt.Tables()
cso_2001 = tbls["2001 CSO Super Preferred Select and Ultimate - Male Nonsmoker, ANB"]
cso_1980 = tbls["1980 CSO - Male Nonsmoker, ANB"]

age = 27
plot([λ -> mt.qx(cso_2001,age,λ),λ -> mt.qx(cso_1980,age,λ)], 1:1:(100-age),label = ["2001 CSO M SuperPref NS" "1980 CSO M NS"], plot_title = ["Comparison of 1980 and 2001 CSO"])
```
![plot of q's](https://i.imgur.com/BvsplkB.png)



### Some Batteries Included

Comes with some tables built in via [mort.SOA.org](https://mort.soa.org) and by using [you agree to their terms](https://mort.soa.org/TermsOfUse.aspx). 

Not all tables have been tested that they work by default, though I have not encountered issues with any of the the VBT/CSO/other usual tables.

Included:
```
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


## License
[![FOSSA Status](https://app.fossa.io/api/projects/git%2Bgithub.com%2Falecloudenback%2FMortalityTables.jl.svg?type=large)](https://app.fossa.io/projects/git%2Bgithub.com%2Falecloudenback%2FMortalityTables.jl?ref=badge_large)