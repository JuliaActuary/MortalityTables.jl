"""
    UltimateMortality(vector; start_age=0)

Given a vector of rates, returns an `OffsetArray` that is indexed by attained age.

Any `AbstractVector` is accepted (a `Vector`, a range, a `view`, or a vector containing `missing`); the input is wrapped without copying.

Give the optional keyword argument to start the indexing at an age other than zero.

# Examples
```julia-repl
julia> m = UltimateMortality([0.1,0.3,0.6,1]);

julia> m[0]
0.1

julia> m = UltimateMortality([0.1,0.3,0.6,1], start_age = 18);

julia> m[18]
0.1

```
"""
function UltimateMortality(v::AbstractVector; start_age = 0)
    # the offset is relative to `v`'s own axes, which need not start at 1 (a view, or a vector
    # that is already indexed by age)
    return OffsetArray(v, start_age - firstindex(v))
end

"""
    SelectMortality(select, ultimate; start_age=0)

Given a matrix of rates, where each row holds the select rates for one issue age, creates an `OffsetArray` that is indexed by issue age, containing a vector of rates indexed by attained age. The ultimate mortality vector is used for rates in the post-select period.

Give the optional keyword argument to start the indexing at an age other than zero.

# Examples
``` 
julia> ult = UltimateMortality([x / 100 for x in 0:100]);

julia> matrix = rand(50,10); # represents random(!) mortality rates with a select period of 10 years

julia> sel = SelectMortality(matrix,ult,start_age=0);

julia> sel[0] # the mortality vector for a select life with issue age 0
 0.12858960119349439
 0.1172480189376135
 0.8237661916705163
 ⋮
 0.98
 0.99
 1.0

julia> sel[0][95] # the mortality rate for a life age 95, that was issued at age 0
0.95
```
"""
function SelectMortality(select, ultimate; start_age = 0)
    # iterate down the rows (issue ages)
    vs = map(enumerate(eachrow(select))) do (i, row)
        _select_row(start_age + i - 1, row, ultimate)
    end

    return OffsetArray(vs, start_age - 1)
end

# One row of a select table: `select_rates` covers attained ages
# `issue_age:issue_age+length(select_rates)-1` (leading `missing` allowed);
# the ultimate table supplies every age after that, through its omega. A
# select period that runs past the ultimate omega simply has no ultimate tail.
function _select_row(issue_age::Integer, select_rates::AbstractVector, ultimate::AbstractVector)
    last_select_age = issue_age + length(select_rates) - 1
    return OffsetArray([select_rates; ultimate[last_select_age+1:end]], issue_age - 1)
end



"""
    MortalityTable(ultimate)
    MortalityTable(select, ultimate)
    MortalityTable(select, ultimate; metadata::TableMetaData)

Constructs a container object which can hold either:
- ultimate-only rates (an `UltimateTable`)
- select and ultimate rates (a `SelectUltimateTable`)

Also pass a keyword argument `metadata=TableMetaData(...)` to store relevant information (source, notes, etc) about the table itself.

# Examples
```julia
# first construct the underlying data
ult = UltimateMortality([x / 100 for x in 0:100]);
matrix = rand(50,10); # random(!) rates for issue ages 0 to 49, with a select period of 10 years
sel = SelectMortality(matrix,ult,start_age=0);

table = MortalityTable(sel,ult)

# can now get rates, indexed by attained age:

table.select[10] # the vector of rates for a risk issued select at age 10 

table.ultimate[99] # 0.99

```
"""
abstract type MortalityTable end

struct SelectUltimateTable{S,U} <: MortalityTable
    select::S
    ultimate::U
    metadata::TableMetaData
end

struct UltimateTable{U} <: MortalityTable
    ultimate::U
    metadata::TableMetaData
end

Base.:(==)(tbl1::UltimateTable, tbl2::UltimateTable) = tbl1.metadata == tbl2.metadata && isequal(tbl1.ultimate, tbl2.ultimate)
function Base.:(==)(tbl1::SelectUltimateTable, tbl2::SelectUltimateTable)
    return (
        tbl1.metadata == tbl2.metadata &&
        isequal(tbl1.ultimate, tbl2.ultimate) &&
        isequal(tbl1.select, tbl2.select)
    )
end
# equal tables (`==`, and so `isequal`) hash alike: the same fields, hashed as `isequal` compares them
Base.hash(tbl::UltimateTable, h::UInt) = hash(tbl.ultimate, hash(tbl.metadata, hash(UltimateTable, h)))
function Base.hash(tbl::SelectUltimateTable, h::UInt)
    return hash(tbl.select, hash(tbl.ultimate, hash(tbl.metadata, hash(SelectUltimateTable, h))))
end



function MortalityTable(select, ultimate; metadata = TableMetaData())
    return SelectUltimateTable(select, ultimate, metadata)
end

function MortalityTable(ultimate; metadata = TableMetaData())
    return UltimateTable(ultimate, metadata)
end


function Base.show(io::IO, ::MIME"text/plain", mt::MortalityTable)
    print(
        io,
        """
        MortalityTable ($(mt.metadata.content_type)):
           Name:
               $(mt.metadata.name)
           Fields:
               $(fieldnames(typeof(mt)))
           Provider:
               $(mt.metadata.provider)
        """,
    )
    # only mort.SOA.org sourced tables have an id and a link
    if mt.metadata.id !== nothing
        print(
            io,
            """
               mort.SOA.org ID:
                   $(mt.metadata.id)
               mort.SOA.org link:
                   https://mort.soa.org/ViewTable.aspx?&TableIdentity=$(mt.metadata.id)
            """,
        )
    end
    print(
        io,
        """
           Description:
               $(mt.metadata.description)
        """,
    )
end


"""
    survival(mortality_vector,to_age)
    survival(mortality_vector,from_age,to_age)

Returns the survival through attained age `to_age`. The start of the calculation is either the start of the vector, or attained_age `from_age`. `from_age` and `to_age` need to be Integers. Add a DeathDistribution as the last argument to handle floating point and non-whole ages:

    survival(mortality_vector,to_age,::DeathDistribution)
    survival(mortality_vector,from_age,to_age,::DeathDistribution)

Survival from a fractional `from_age` is conditional on surviving to `from_age`: it equals `survival(v, to_age, dd) / survival(v, from_age, dd)` under the same assumption, so survival over consecutive intervals multiplies. Where the assumption gives zero survival to a fractional `from_age` (after a rate of one under `ConstantForce` or `Balducci`) there is nothing to condition on; the formulas are still evaluated as written and give finite values.

When `to_age` is before `from_age`, the result is the reverse factor `1 / survival(v, to_age, from_age)`: the number expected alive at the earlier age for each life alive at the later one, as used to project a population backward or to accumulate with the benefit of survivorship. It is not a probability (it can exceed one, and the corresponding `decrement` is negative), and it is an expected-value back-calculation rather than a reconstruction of realized deaths. It is defined where the forward survival is positive; a zero forward survival gives `Inf`. With it, survival composes over any three ages whose factors are positive and representable: `survival(v, a, c) == survival(v, a, b) * survival(v, b, c)`, whatever their order. (A zero factor has no inverse: `0 * Inf` is not one.) Ages outside the table (such as a negative `to_age` for a table starting at zero) are a `BoundsError`.

Results have the numeric type of the rates, promoted with the ages' type for fractional ages, including the exact one of an empty interval: a `BigFloat` table gives `BigFloat` survival, and a `Float32` table gives `Float32` survival at whole ages. A vector whose element type doesn't name the rates' number type (such as `Any`, `Real`, or `Missing` for a column with no rates) has `Float64` identities, as in version 2: an empty interval gives `1.0`, and a nonempty result takes the type its rates' arithmetic gives. Use a concretely typed vector for type-consistent and fast results.

Parametric models (see `ParametricMortality`) are continuous and need no fractional-age assumption, so they accept a trailing `DeathDistribution` and ignore it.

# Examples
```julia-repl
julia> qs = UltimateMortality([0.1,0.3,0.6,1]);
    
julia> survival(qs,0)
1.0
julia> survival(qs,1)
0.9

julia> survival(qs,1,1)
1.0
julia> survival(qs,1,2)
0.7

julia> survival(qs,0.5,UniformDeaths())
0.95
```
"""
function survival(v::AbstractArray, to_age)
    return survival(v, firstindex(v), to_age)
end
function survival(v::AbstractArray, to_age, dd::DeathDistribution)
    return survival(v, firstindex(v), to_age, dd)
end

_decrement(surv, q) = surv * (1 - q)

# The numeric type of survival and decrement: that of the rates (ignoring `missing`, as in a
# select row with a gap), promoted with the ages' type for fractional ages. Identities such as
# the survival over an empty interval have this type too, so a BigFloat or Float32 table keeps
# its precision whether or not an interval is empty. An element type that isn't concrete says
# nothing about the rates' numbers (`Any`, `Real`, or `Union{}` for an all-`missing` vector), so
# the identities are Float64's, as before v3. A concrete non-number, such as a `String` rate,
# still fails where it is converted.
_rate_type(v) = _rate_type(nonmissingtype(eltype(v)))
_rate_type(::Type{T}) where {T} = isconcretetype(T) ? float(T) : Float64
_survival_type(v, ages...) = float(promote_type(_rate_type(v), map(typeof, ages)...))

function survival(v::AbstractArray, from_age::Int, to_age::Int)
    # a reversed interval is the reverse factor (see the docstring)
    from_age > to_age && return inv(survival(v, to_age, from_age))
    # an empty age range (from_age == to_age) reduces to `init`, i.e. one
    return @views reduce(_decrement, v[from_age:(to_age-1)], init = one(_rate_type(v)))
end

function survival(v::AbstractArray, from_age, to_age, dd::DeathDistribution)
    # a reversed interval is the reverse factor (see the docstring)
    from_age > to_age && return inv(survival(v, to_age, from_age, dd))
    T = _survival_type(v, from_age, to_age)
    # the survival over the whole ages, times the partial years before and after them
    age_low = ceil(Int, from_age)
    age_high = floor(Int, to_age)
    # within one year of age
    age_high < age_low && return 1 - decrement_partial_year(v, from_age, to_age, dd)
    low_residual = age_low == from_age ? one(T) : 1 - decrement_partial_year(v, from_age, age_low, dd)
    high_residual = age_high == to_age ? one(T) : 1 - decrement_partial_year(v, age_high, to_age, dd)
    # an empty interval reduces to `init`, i.e. one
    whole = @views reduce(_decrement, v[age_low:(age_high-1)], init = one(T))
    return whole * low_residual * high_residual
end

# Reference: Experience Study Calculations, 2016, Society of Actuaries
# https://www.soa.org/globalassets/assets/Files/Research/2016-10-experience-study-calculations.pdf
#
# The decrement between `from_age` and `to_age` within one year of age x = ⌊from_age⌋,
# conditional on surviving to `from_age`. With s = from_age - x and t = to_age - x, it is
# 1 - S(x+t)/S(x+s), where S(x+u)/S(x) is 1 - u·q (UniformDeaths), (1-q)^u (ConstantForce), or
# (1-q)/(1-(1-u)·q) (Balducci).
function decrement_partial_year(v, from_age, to_age, dd::UniformDeaths)
    x = floor(Int, from_age)
    q = v[x]
    s, t = from_age - x, to_age - x
    return q * (t - s) / (1 - s * q)
end

function decrement_partial_year(v, from_age, to_age, dd::ConstantForce)
    return 1 - (1 - v[floor(Int, from_age)])^(to_age - from_age)
end

function decrement_partial_year(v, from_age, to_age, dd::Balducci)
    x = floor(Int, from_age)
    q = v[x]
    s, t = from_age - x, to_age - x
    return q * (t - s) / (1 - q + t * q)
end

"""
    decrement(mortality_vector,to_age)
    decrement(mortality_vector,from_age,to_age)

Returns the cumulative decrement through attained age `to_age`. The start of the calculation is either the start of the vector, or attained_age `from_age`. `from_age` and `to_age` need to be Integers. Add a DeathDistribution as the last argument to handle floating point and non-whole ages:

    decrement(mortality_vector,to_age,::DeathDistribution)
    decrement(mortality_vector,from_age,to_age,::DeathDistribution)

# Examples
```julia-repl
julia> qs = UltimateMortality([0.1,0.3,0.6,1]);
    
julia> decrement(qs,0)
0.0
julia> decrement(qs,1)
0.1

julia> decrement(qs,1,1)
0.0
julia> decrement(qs,1,2)
0.3

julia> decrement(qs,0.5,UniformDeaths())
0.05
```

The decrement is accumulated directly (``d \\leftarrow q + d(1 - q)``) rather than computed as one minus the survival product, so a small decrement keeps its precision. A reversed interval gives the negative decrement `1 - survival(v, from_age, to_age)` of the reverse factor, computed as `-d / S` from the forward decrement `d` and forward survival `S`, so it keeps its precision both for a small decrement and for a small forward survival.

A type that defines only `survival` gets `decrement` as `1 - survival`; define `decrement` too where that complement would lose precision.
"""
decrement(v::AbstractArray, to_age) = decrement(v, firstindex(v), to_age)
decrement(v::AbstractArray, to_age, dd::DeathDistribution) = decrement(v, firstindex(v), to_age, dd)

# the extension contract: a type that defines `survival` has the complementary `decrement`
decrement(v, args...) = 1 - survival(v, args...)

# d ← q + d·(1 - q), i.e. 1 - (1 - d)(1 - q): the complement of the survival product without
# the cancellation in 1 - ∏(1 - q), so a small decrement is not rounded away. For rates in
# [0, 1] both terms are non-negative, so the sum is accurate to a few ulps; `1 - q` does not
# depend on d, so each step is one fused multiply-add, as fast as the product. Values above one
# (claim costs or factors, which some bundled tables hold) are not probabilities; the same
# expression is evaluated for them, and it can round more than the product would.
_accumulate_decrement(d, q) = muladd(d, 1 - q, q)

# The decrement of a reversed interval, 1 - 1/S = -d/S for the forward decrement d and survival
# S. Both are accurate to a few ulps, so the quotient is too: rebuilding S as 1 - d would round it
# to zero once d rounds to one, although S (say 2^-60) is representable.
_reverse_decrement(v, args...) = -decrement(v, args...) / survival(v, args...)

function decrement(v::AbstractArray, from_age::Int, to_age::Int)
    from_age > to_age && return _reverse_decrement(v, to_age, from_age)
    return @views reduce(_accumulate_decrement, v[from_age:(to_age-1)], init = zero(_rate_type(v)))
end

function decrement(v::AbstractArray, from_age, to_age, dd::DeathDistribution)
    from_age > to_age && return _reverse_decrement(v, to_age, from_age, dd)
    T = _survival_type(v, from_age, to_age)
    from_age == to_age && return zero(T)
    age_low = ceil(Int, from_age)
    age_high = floor(Int, to_age)
    # within one year of age
    age_high < age_low && return _decrement_piece(v, from_age, to_age, dd)
    d = age_low == from_age ? zero(T) : _decrement_piece(v, from_age, age_low, dd)
    d = @views reduce(_accumulate_decrement, v[age_low:(age_high-1)], init = d)
    return age_high == to_age ? d : _accumulate_decrement(d, _decrement_piece(v, age_high, to_age, dd))
end

# The decrement over part of one year of age: `decrement_partial_year`, except that the
# constant force is written as -expm1((t - s)·log1p(-q)), since 1 - (1 - q)^(t - s) rounds a
# small rate away. (`survival` keeps `1 - decrement_partial_year`, so its values are unchanged.)
_decrement_piece(v, from_age, to_age, dd::DeathDistribution) = decrement_partial_year(v, from_age, to_age, dd)
_decrement_piece(v, from_age, to_age, ::ConstantForce) = -expm1((to_age - from_age) * log1p(-v[floor(Int, from_age)]))

"""
    omega(x)
    ω(x)

Returns the last index of the given vector. For mortality vectors this means the last attained age for which a rate is defined.

Note that `omega` can vary depending on the issue age for a select table, and that a select `omega` may differ from the table's ultimate `omega`.

For a parametric model (see `ParametricMortality`), `omega` is the last age at which its law is defined: `Inf` for a law defined at every age, `m` for `Wittstein` (where its formula ends), and `n` for `VanderMaen` and `VanderMaen2` (whose hazard has a pole there). An age past `omega` is a `DomainError`. Survival need not reach zero at `omega` (for `Wittstein` it does not), so a projection that stops at `omega` truncates such a law. `omega` does not validate a law's parameters either: some parameters give a negative hazard at ages inside the domain.

ω is aliased to omega, but un-exported. To use, do `using MortalityTables: ω` when importing or call `MortalityTables.ω()`

# Examples

```julia-repl
julia> qs = UltimateMortality([0.1,0.3,0.6,1]);
julia> omega(qs)
3

julia> qs = UltimateMortality([0.1,0.3,0.6,1],start_age=10);
julia> omega(qs)
13

```
"""
function omega(x)
    return lastindex(x)
end

const ω = omega


"""
    mortality_vector(vec; start_age=0)

An alias for [`UltimateMortality`](@ref): wraps `vec` in an `OffsetArray` indexed by attained age, starting at `start_age`, rather than always starting from `1`. The package and JuliaActuary ecosystem assume that the rates are indexed by attained age, and this allows transformation of tables without a direct dependency on **OffsetArrays.jl**.
"""
const mortality_vector = UltimateMortality