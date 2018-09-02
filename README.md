# MortalityTables

#### Code Review: [![Build Status](https://travis-ci.org/alecloudenback/MortalityTables.jl.svg?branch=master)](https://travis-ci.org/alecloudenback/MortalityTables.jl) [![Coverage Status](https://coveralls.io/repos/github/alecloudenback/MortalityTables.jl/badge.svg?branch=master)](https://coveralls.io/github/alecloudenback/MortalityTables.jl?branch=master) [![codecov.io](http://codecov.io/github/alecloudenback/MortalityTables.jl/coverage.svg?branch=master)](http://codecov.io/github/alecloudenback/MortalityTables.jl?branch=master)
A Julia package for working with MortalityTables. Has first-class support for missing values.

Currently Available: `2001 CSO` and `2001 VBT`

### Usage Example

```julia
import MortalityTables as mt

tables = mt.Tables()
cso = tables["2001 CSO Super Preferred Select and Ultimate - Male Nonsmoker, ANB"]


mt.qx(cso,35,1) # == .00037
mt.qx(cso,35,61) # == .26719
mt.qx(cso,95) # == .26719
mt.qx(cso,35,95) # == missing (table doesn't extend that far)
```


Comes with some tables built in.

Todo:
- Docs
- Automatically parse built-in tables
- Add more built-in tables
- Usage Examples
- More tests
- Performance testing
