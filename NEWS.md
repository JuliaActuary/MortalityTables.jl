# MortalityTables.jl release notes

## v3.0.0

### Breaking changes

- Julia 1.10 or newer is required.
- Indexing a parametric model (`m[x]`) is removed. Use `hazard(m, x)` or `m(x)`.
- `UltimateTable` no longer forwards indexing or `lastindex`; index the rates directly with `table.ultimate[age]`. Tables are opaque containers: passing a whole table to `survival`, `decrement`, or `life_expectancy` is a `MethodError` (2.x threw an `ArgumentError` pointing to `table.ultimate`).
- `survival` and `decrement` on rate vectors now require an `AbstractArray` first argument. This removes the method ambiguities that previously turned some misuse into an ambiguity error.
- A reversed interval of a rate vector (`to_age` before `from_age`) is the reverse factor `1 / survival(v, to_age, from_age)`, as it already was for parametric laws and is for FinanceCore and FinanceModels discount factors. 2.x returned `1.0` for whole ages and a meaningless value for fractional ones. The reverse factor projects a population backward or accumulates with the benefit of survivorship; it is not a probability (it can exceed one, and its `decrement` is negative), and survival now composes over any three ages whose factors are positive. Ages outside a table are a `BoundsError`: a negative `to_age` no longer returns `1.0`.
- `survival` and `decrement` on a rate vector return the rates' numeric type, promoted with the ages' type for fractional ages, including the exact one and zero of an empty interval. A `Float32` table gives `Float32` results at whole ages (2.x widened them to `Float64`), and a `BigFloat` table gives `BigFloat` results for empty intervals too (2.x returned the `Float64` literal `1.0`). `Float64` tables are unchanged, and so is a vector whose element type isn't concrete (`Any`, `Real`, or `Missing` for a column with no rates): its identities stay the `Float64` `1.0` and `0.0`, and a nonempty result has the type of its rates' arithmetic.
- Parametric laws have concrete, promoted field types (`Makeham{Float64}`). Constructors accept any mix of `Real` keywords (`Makeham(a=1, b=0.5, c=0)` is a `Makeham{Float64}`), and every default law is `isbits`. `survival` is roughly 20x faster for laws without a closed-form cumulative hazard.
- `MakehamGompertz` was removed from the export list (it was never defined).
- Loading a table from a CSV file is a package extension activated by `using CSV`; the call `MortalityTable(CSV.File(...))` is unchanged.

### Bug fixes

- `survival` and `decrement` on a rate vector from a fractional starting age are now conditional on surviving to that age under `Uniform()` and `Balducci()`. In 2.x the partial year used only the interval's length: with `q = 0.2`, survival from age 0.5 to 1 was `0.9` under `Uniform()` instead of `0.8/0.9`, and two half years multiplied to `0.81` instead of the year's `0.8`. Values change for intervals that start between birthdays; intervals that start on a birthday, and `Constant()`, are unchanged.
- Reading a CSV table places its rates by their age and duration labels, as XTbML now does. A table with grouped ages, gaps, or nonconsecutive duration columns no longer shifts rates to other ages (rates labelled ages 60 and 62 loaded as ages 60 and 61). A select row with no rates is an absent issue age (`missing`), as in XTbML, instead of a zero-length select period that used the ultimate rates.
- `Beard`, `MakehamBeard`, `Kannisto`, `KannistoMakeham`, `GammaGompertz`, and `Martinelle` no longer return `NaN` where `exp(b·age)` overflows; their bounded hazards stay finite. For example `survival(Kannisto(a = 0.5, b = 10), 80, 81)` is `exp(-1)` rather than `NaN`. The algebra follows the sign of the exponent, so hazards in the negative-growth tail keep their values too (including subnormal ones, and `GammaGompertz` with zero frailty). Derivatives with respect to each parameter stay finite as well, including at a zero amplitude, where the stable forms keep the laws' smooth `exp(b·age)` sensitivity.
- Kannisto's `cumhazard`, and therefore its `survival`, returned wrong values in 2.x (`survival(Kannisto(), 0)` was 1.76). The closed form had its `a` and `b` parameters swapped.
- Two-age `survival` on a parametric model is now the difference of cumulative hazards rather than a ratio of survival probabilities, so it is well defined where survival underflows to zero (previously `NaN`).
- XTbML rates are placed by their age and duration labels instead of by position. Eight mort.SOA.org tables label grouped ages; for example, t352 ("1946-49 Basic Table, ANB") has issue ages 12, 17, 22, …. For these tables most ages returned another age's rates: `select[17]` gave issue age 37's. Each rate now sits at its own age, with `missing` between the groups. An empty cell inside a select row no longer shifts the later rates to earlier ages, and a repeated label throws.
- `Makeham`, `Gompertz`, and `Kannisto` with `b = 0` return the constant-hazard limit instead of `NaN`. At an infinite age a zero coefficient no longer turns the result into `NaN`: `survival(Gompertz(), Inf)` and a constant hazard's survival there are `0`, a law with no hazard keeps survival `1`, and `hazard(Makeham(b = 0), Inf)` is `a + c`.
- 34 additional mort.SOA.org tables with an empty metadata element now load (for example `table(217)`). The 760 bundled tables that still do not load (other layouts, such as claim termination rates or projection scales, and select durations labelled from 0) are listed in `test/data/unsupported_tables.txt`.
- An XTbML or CSV table with more than two tables of rates raises an `ArgumentError` naming the table, rather than loading with its middle table dropped. Four bundled tables have three (t357, t359, t754, and t755); they did not load in 2.x either.
- The `get_SOA_table` "not found" error message now includes the requested table name.
- Showing a custom table no longer prints a broken mort.SOA.org link.

### Other changes

- `cumhazard` is defined for every parametric law. It integrates `hazard` numerically unless the law provides a closed form, and `survival` follows from it.
- Parametric models accept and ignore a `DeathDistribution` in `survival`, `decrement`, and `life_expectancy`, so a model works wherever a table's rate vector does. `omega` of a parametric model is the last age at which its law is defined: `Inf` for most laws, `m` for `Wittstein`, and `n` for `VanderMaen` and `VanderMaen2`. Survival need not be zero there, so a projection over a parametric law needs its own horizon.
- `decrement` on a parametric model is `-expm1(-cumhazard(...))`, and on a rate vector the decrement is accumulated directly (`d ← d + q(1 - d)`) rather than computed as one minus the survival product, so a small decrement is not rounded away: a rate of `1e-18` gave `0.0`, and `1e-6` kept only about 11 digits. Survival values are unchanged. A type that defines only `survival` still gets `decrement` as `1 - survival`.
- `life_expectancy` on a rate vector is linear rather than quadratic in the number of remaining ages and allocation free. An age past omega is a `BoundsError`.
- `UltimateMortality` accepts any `AbstractVector` (ranges, views, vectors containing `missing`, or a vector already indexed by age) and wraps it without copying; `start_age` sets the first age whatever the input's own axes. `mortality_vector` is an alias for it.
- XTbML parsing uses [XML.jl](https://github.com/JuliaData/XML.jl): loading a table is about 2 to 3 times faster and `using MortalityTables` no longer loads libxml2.
- `Pkg`, `UnPack`, `Requires`, `Memoize`, `Parsers`, and `XMLDict` are no longer dependencies.
- The Dukes-MacDonald functions are documented on the Tables page.
