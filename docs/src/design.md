# The design of MortalityTables.jl

## MortalityTables.jl types

There are two types of MortalityTable subtypes for SOA data: `SelectUltimateTable` and `UltimateTable`, so define a function for each. 

A table with both select and ultimate rates:

```@example
MortalityTables.table("2001 VBT Residual Standard Select and Ultimate - Male Nonsmoker, ANB")
```

A table with just ultimate rates:
```@example
MortalityTables.table("1941 CSO Basic Table, ANB")
```

## Motivation for the current design
x
### Encourage data-driven design

Probably best explained with an example. In the following example, we calculate the ultimate survivorship for every age in the `cso` table from above. The table itself gives us the right range of ages to do this with. No need to assume a 100 or 121 omega age and handle edge cases when working with tables that don't line up to the usual expectations.

```julia
map(eachindex(cso.ultimate)) do issue_age
	survival(cso.ultimate,issue_age,omega(cso.ultimate))
end
```

Of course, here we are limited by the quality of the input tables. For example, in the `vbt` table from the SOA loaded above, select rates are missing from under-18s, even though select rates are defined for a life issued under-18.

### Speed of calculation

In modeling contexts, efficiency is often a key part of the model. This is probably close to the fastest single-threaded computation of survival probabilities possible, and easy to extend to parallel calculations.

### Parametric models

With MortalityTables.jl, virtually transparent to the user of the package, you can drop in parametric models in place of tabular-SOA-like tables. Here is a Makeham mortality model and the same syntax is used as with the SOA tables:

```julia
survival(MortalityTables.Makeham(),25,50) # 25 to 50 year old survival
```

### Table metadata

When you display a table, you get to see related metadata, which can be lost if simply parsing into a matrix or dataframe:

```@example
MortalityTables.table("2001 VBT Residual Standard Select and Ultimate - Male Nonsmoker, ANB")
```

### Partial Year assumptions

Built in are partial year assumptions, making it easy to use more realistic assumptions like Uniform death distribution. The types `Balducci()`, `Uniform()`, and `Constant()` are used to dispatch to the appropriate function without any runtime overhead.

## Where MortalityTables.jl falls short

In no particular order:

- it mainly relies on data from mort.SOA.org. Some of the tables have issues like inconsistent naming, missing fields
- would be nice to offer 'sets' of table to users. Early tests failed to do this in an automated way because of issues with the naming consistency of tables mentioned in the first bullet
- (currently) no automatic conversion to a dataframe-compatible, though this is discussed above
- MortalityTables.jl is a misnomer, as the `ORT` of MORT.SOA.org is 'other rate tables' which are offered in MortalityTables.jl
- no built-in support for mortality improvement
- has issues with autodifferentiation because of the discrete and non-continuous nature of the annual rate tables