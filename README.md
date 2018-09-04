# MortalityTables

#### Code Review: [![Build Status](https://travis-ci.org/alecloudenback/MortalityTables.jl.svg?branch=master)](https://travis-ci.org/alecloudenback/MortalityTables.jl) [![Coverage Status](https://coveralls.io/repos/github/alecloudenback/MortalityTables.jl/badge.svg?branch=master)](https://coveralls.io/github/alecloudenback/MortalityTables.jl?branch=master) [![codecov.io](http://codecov.io/github/alecloudenback/MortalityTables.jl/coverage.svg?branch=master)](http://codecov.io/github/alecloudenback/MortalityTables.jl?branch=master)
A Julia package for working with MortalityTables. Has first-class support for missing values.

Currently Available: `2001 CSO` and `2001 VBT`

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


Comes with some tables built in. Currently these are included:

```

2001 CSO Super Preferred Select and Ultimate - Male Nonsmoker, ANB
2001 CSO Preferred Select and Ultimate - Male Nonsmoker, ANB
2001 CSO Residual Standard Select and Ultimate  - Male Nonsmoker, ANB
2001 CSO Preferred Select and Ultimate - Male Smoker, ANB
2001 CSO Residual Standard Select and Ultimate - Male Smoker, ANB
2001 CSO Super Preferred Select and Ultimate - Female Nonsmoker, ANB
2001 CSO Preferred Select and Ultimate - Female Nonsmoker, ANB
2001 CSO Residual Standard Select and Ultimate - Female Nonsmoker, ANB
2001 CSO Preferred Select and Ultimate - Female Smoker, ANB
2001 CSO Residual Standard Select and Ultimate - Female Smoker, ANB
2001 CSO Super Preferred Select and Ultimate - Male Nonsmoker, ALB
2001 CSO Preferred Select and Ultimate - Male Nonsmoker, ALB
2001 CSO Residual Standard Select and Ultimate - Male Nonsmoker, ALB
2001 CSO Preferred Select and Ultimate - Male Smoker, ALB
2001 CSO Residual Standard Select and Ultimate - Male Smoker, ALB
2001 CSO Super Preferred Select and Ultimate - Female Nonsmoker, ALB
2001 CSO Preferred Select and Ultimate - Female Nonsmoker, ALB
2001 CSO Residual Standard Select and Ultimate - Female Nonsmoker, ALB
2001 CSO Preferred Select and Ultimate - Female Smoker, ALB
2001 CSO Residual Standard Select and Ultimate - Female Smoker, ALB
2001 VBT Super Preferred Select and Ultimate - Male Nonsmoker, ANB
2001 VBT Preferred Select and Ultimate - Male Nonsmoker, ANB
2001 VBT Residual Standard Select and Ultimate - Male Nonsmoker, ANB
2001 VBT Preferred Select and Ultimate - Male Smoker, ANB
2001 VBT Residual Standard Select and Ultimate - Male Smoker, ANB
2001 VBT Super Preferred Select and Ultimate - Female Nonsmoker, ANB
2001 VBT Preferred Select and Ultimate - Female Nonsmoker, ANB
2001 VBT Residual Standard Select and Ultimate - Female Nonsmoker, ANB
2001 VBT Preferred Select and Ultimate - Female Smoker, ANB
2001 VBT Residual Standard Select and Ultimate - Female Smoker, ANB
2001 VBT Super Preferred Select and Ultimate - Male Nonsmoker, ALB
2001 VBT Preferred Select and Ultimate - Male Nonsmoker, ALB
2001 VBT Preferred Select and Ultimate - Male Smoker, ALB
2001 VBT Residual Standard Select and Ultimate - Male Nonsmoker, ALB
2001 VBT Residual Standard Select and Ultimate - Male Smoker, ALB
2001 VBT Super Preferred Select and Ultimate - Female Nonsmoker, ALB
2001 VBT Preferred Select and Ultimate - Female Nonsmoker, ALB
2001 VBT Residual Standard Select and Ultimate - Female Nonsmoker, ALB
2001 VBT Preferred Select and Ultimate - Female Smoker, ALB
2001 VBT Residual Standard Select and Ultimate - Female Smoker, ALB
2001 CSO Select and Ultimate – Male Composite, ANB
2001 CSO Select and Ultimate - Male Nonsmoker, ANB
2001 CSO Select and Ultimate  - Male Smoker, ANB
2001 CSO Select and Ultimate - Female Composite, ANB
2001 CSO Select and Ultimate - Female Nonsmoker, ANB
2001 CSO Select and Ultimate - Female Smoker, ANB
2001 CSO Composite Select and Ultimate - Male, ALB
2001 CSO Composite Select and Ultimate - Female, ALB
2001 CSO Select and Ultimate - Male Nonsmoker, ALB
2001 CSO Select and Ultimate - Female Nonsmoker, ALB
2001 CSO Select and Ultimate  - Male Smoker, ALB
2001 CSO Select and Ultimate - Female Smoker, ALB
1980 CSO Basic Table – Female, ANB
1980 CSO Basic Table - Female Nonsmoker, ANB
1980 CSO Basic Table - Female Smoker, ANB
1980 CSO Basic Table – Male, ANB
1980 CSO Basic Table - Male Nonsmoker, ANB
1980 CSO Basic Table - Male Smoker, ANB
1980 CET – Female, ALB
1980 CET - Female Nonsmoker, ALB
1980 CET - Female Nonsmoker, ANB
1980 CET - Female Smoker, ALB
1980 CET - Female Smoker, ANB
1980 CET – Male, ALB
1980 CET - Male Nonsmoker, ALB
1980 CET - Male Nonsmoker, ANB
1980 CET - Male Smoker, ALB
1980 CET  - Male Smoker, ANB
1980 CSO – Female, ALB
1980 CSO - Female, ANB
1980 CSO - Female Nonsmoker, ALB
1980 CSO - Female Nonsmoker, ANB
1980 CSO - Female Smoker, ALB
1980 CSO - Female Smoker, ANB
1980 CSO – Male, ALB
1980 CSO  - Male, ANB
1980 CSO - Male Nonsmoker, ALB
1980 CSO - Male Nonsmoker, ANB
1980 CSO - Male Smoker, ALB
1980 CSO - Male Smoker, ANB
```

Todo:
- Docs
- Automatically parse built-in tables
- Add more built-in tables
- Usage Examples
- More tests
- Performance testing
