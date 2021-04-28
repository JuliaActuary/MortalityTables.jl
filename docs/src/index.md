# MortalityTables.jl

```@meta
DocTestSetup = quote
    using MortalityTables
end
```

Part of the [JuliaActuary.org](https://juliaactuary.org/) effort to build packages that enable actuaries everywhere to build solutions using open-source tools.

## Overview

A Julia package for working with MortalityTables. Has:

- Lots of bundled SOA mort.soa.org tables
- `survival` and `decrement` functions to calculate decrements over period of time
- Partial year mortality calculations (Uniform, Constant, Balducci)
- Friendly syntax and flexible usage
- Extensive set of parametric mortality models.

## Installation

Installation
The package can be installed with the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia-repl
pkg> add MortalityTables
```

Or, equivalently, via the Pkg API:

```julia-repl
julia> import Pkg; Pkg.add("MortalityTables")
```

## Quickstart

Loading the package and bundled tables:

```julia
julia> using MortalityTables

julia> tables = MortalityTables.tables()
Dict{String,MortalityTable} with 266 entries:
  "2015 VBT Female Non-Smoker RR90 ALB"                                       => SelectUltimateTable{OffsetArray{OffsetArray{Float64,1,Array{Float64,1}},1,Array{OffsetArray{F…  
  "2017 Loaded CSO Preferred Structure Nonsmoker Preferred Female ANB"        => SelectUltimateTable{OffsetArray{OffsetArray{Float64,1,Array{Float64,1}},1,Array{OffsetArray{F…  
  ⋮                                                                            => ⋮
```

Get information about a particular table:

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

Calculate the force of mortality or survival over a range of time:

```julia
julia> survival(vbt2001.ultimate,30,40) # the survival between ages 30 and 40
0.9894404665434904

julia> decrement(vbt2001.ultimate,30,40) # the decrement between ages 30 and 40
0.010559533456509618
```

Non-whole periods of time are supported when you specify the assumption (`Constant()`, `Uniform()`, or `Balducci()`) for fractional periods:

```julia
julia> survival(vbt2001.ultimate,30,40.5,Uniform()) # the survival between ages 30 and 40.5
0.9887676470262408
```

## License and Usage Notes

This package is [MIT Licenced](https://github.com/JuliaActuary/MortalityTables.jl/blob/master/LICENSE) ([tl;dr explanation](https://tldrlegal.com/license/mit-license)).

This package contains tables downloaded from [mort.SOA.org](https://mort.soa.org). By using those tables within this package, you also agree to [their Terms of Use](https://mort.soa.org/TermsOfUse.aspx).

## Discussion and Questions

If you have other ideas or questions, feel free to also open an issue, or discuss on the community [Zulip](https://julialang.zulipchat.com/#narrow/stream/249536-actuary) or [Slack #actuary channel](https://julialang.org/slack/). We welcome all actuarial and related disciplines!
