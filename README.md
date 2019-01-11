# MortalityTables

#### Code Review: [![Build Status](https://travis-ci.org/alecloudenback/MortalityTables.jl.svg?branch=master)](https://travis-ci.org/alecloudenback/MortalityTables.jl) [![Coverage Status](https://coveralls.io/repos/github/alecloudenback/MortalityTables.jl/badge.svg?branch=master)](https://coveralls.io/github/alecloudenback/MortalityTables.jl?branch=master) [![codecov.io](http://codecov.io/github/alecloudenback/MortalityTables.jl/coverage.svg?branch=master)](http://codecov.io/github/alecloudenback/MortalityTables.jl?branch=master)
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


[Click here to see which tables are included.](BundledTables.md)



### Todos

- Docs
- Automatically parse built-in tables
- Add more built-in tables
- Usage Examples
- More tests
- Performance testing
