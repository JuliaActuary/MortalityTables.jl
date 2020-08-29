"""
    UltimateMortality(vector; start_age=0)

Given a vector of rates, returns an `OffsetArray` that is indexed by attained age. 

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
function UltimateMortality(v::Array{<:Real,1}; start_age = 0)
    return OffsetArray(v, start_age - 1)
end

"""
    SelectMortality(select, ultimate; start_age=0)

Given a matrix rates, where the first row represents the select rates for a risk, will create a an `OffsetArray` that is indexed by issue age, containing a vector of rate indexed by attained age. The ultimate mortality vector is used for rates in the post-select period.

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
    vs = map(enumerate(eachrow(select))) do (row_num, row)
        end_age = start_age + (row_num - 1) + (length(row) - 1)
        OffsetArray([row ; ultimate[end_age + 1:end]], (start_age - 1) + (row_num - 1))
    end

    return OffsetArray(vs, start_age - 1)
end



"""
    MortalityTable(ultimate)
    MortalityTable(select, ultimate)
    MortalityTable(select, ultimate; metadata::MetaData)

Constructs a container object which can hold either:
    - ultimate-only rates (an `UltimateTable`)
    - select and ultimate rates (a `SelectUltimateTable`)

Also pass a keyword argument `metadata=MetaData(...)` to store relevant information (source, notes, etc) about the table itself.

# Examples
```julia
# first construct the underlying data
ult = UltimateMortality([x / 100 for x in 0:100]); # first ma
matrix = rand(10,50); # represents random mortality rates with a select period of 10 years
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

Base.getindex(u::UltimateTable,x) = u.ultimate[x]
Base.lastindex(u::UltimateTable) = lastindex(u.ultimate)


function MortalityTable(select, ultimate; metadata = TableMetaData())
    return SelectUltimateTable(select, ultimate, metadata)
end

function MortalityTable(ultimate; metadata = TableMetaData())
    return UltimateTable(ultimate, metadata)
end


Base.show(io::IO, ::MIME"text/plain", mt::MortalityTable) = print(
    io,
    """
    MortalityTable ($(mt.metadata.content_type)):
       Name:
           $(mt.metadata.name)
       Fields: 
           $(fieldnames(typeof(mt)))
       Provider:
           $(mt.metadata.provider)
       mort.SOA.org ID:
           $(mt.metadata.id)
       mort.SOA.org link:
           https://mort.soa.org/ViewTable.aspx?&TableIdentity=$(mt.metadata.id)
       Description:
           $(mt.metadata.description)
    """,
)


"""
    survivorship(mortality_vector,to_age)
    survivorship(mortality_vector,from_age,to_age)

Returns the survivorship through attained age `to_age`. The start of the calculation is either the start of the vector, or attained_age `from_age`. `from_age` and `to_age` need to be Integers. Add a DeathDistribution as the last argument to handle floating point and non-whole ages:

    survivorship(mortality_vector,to_age,::DeathDistribution)
    survivorship(mortality_vector,from_age,to_age,::DeathDistribution)

If given a negative `to_age`, it will return `1.0`. Aside from simplifying the code, this makes sense as for something to exist in order to decrement in the first place, it must have existed and surived to the point of  being able to be decremented.

# Examples
```julia-repl
julia> qs = UltimateMortality([0.1,0.3,0.6,1]);
    
julia> survivorship(qs,0)
1.0
julia> survivorship(qs,1)
0.9

julia> survivorship(qs,1,1)
1.0
julia> survivorship(qs,1,2)
0.7

julia> survivorship(qs,0.5,Uniform())
0.95
```
"""
function survivorship(v, to_age)
    return survivorship(v, firstindex(v), to_age)
end
function survivorship(v, to_age, dd::DeathDistribution)
    return survivorship(v, firstindex(v), to_age, dd)
end

function survivorship(v::T, from_age::Int, to_age::Int) where {T <: AbstractArray}
    if from_age == to_age
        return 1.0
    else
        return @views reduce(*,
            1 .- v[from_age:(to_age - 1)],
            init = 1.0
            )
    end
end

function survivorship(v::T, from_age, to_age, dd::DeathDistribution) where {T <: AbstractArray}
    # calculate the survivorship for the rounded ages, and then the high and low high_residual
    age_low = ceil(Int, from_age)
    age_high = floor(Int, to_age)

    #if from_age and to_age are fractional parts of the same attained age, then age_high will round down to 
    # be below the rounded-up age_low. This line will short circuit the rest and just return the fractional year survivorship
    age_high < age_low && return 1 - decrement_partial_year(v, from_age, to_age, dd)
    
    
    if age_low == from_age
        low_residual = 1.0
    else
        low_residual = 1 - decrement_partial_year(v, from_age, age_low, dd)
    end

    if age_high == to_age
        high_residual = 1.0
    else
        high_residual = 1 - decrement_partial_year(v, age_high, to_age, dd)
    end

    if from_age == to_age
        return 1.0
    else
         
        whole = @views reduce(*,
            1 .- v[age_low:(age_high - 1)],
            init = 1.0
            )

        return whole * low_residual * high_residual
    end
end

# Reference: Experience Study Calculations, 2016, Society of Actuaries
# https://www.soa.org/globalassets/assets/Files/Research/2016-10-experience-study-calculations.pdf
function decrement_partial_year(v, from_age, to_age, dd::Uniform)
    return v[floor(Int,from_age)] * (to_age - from_age)
end

function decrement_partial_year(v, from_age, to_age, dd::Constant)
    return 1 - (1 - v[floor(Int,from_age)])^(to_age - from_age)
end

function decrement_partial_year(v, from_age, to_age, dd::Balducci)
    q′ = v[floor(Int,from_age)]
    frac = (to_age - from_age)
    return 1 - (1 - q′) / (1 - (1 - frac) * q′)
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

julia> decrement(qs,0.5,Uniform())
0.05
```
"""
decrement(v,to_age) = 1 .- survivorship(v, to_age)
decrement(v,to_age,dd::DeathDistribution) = 1 .- survivorship(v, to_age, dd)  
decrement(v,from_age,to_age) = 1 .- survivorship(v, from_age, to_age) 
decrement(v,from_age,to_age,dd::DeathDistribution) = 1 .- survivorship(v, from_age, to_age, dd) 


"""
    omega(x)
    ω(x)

Returns the last index of the given vector. For mortality vectors this means the last attained age for which a rate is defined.

ω is aliased to omega, but unexported. To use, do `using MortalityTables: ω` when importing or call `MortalityTables.ω()`

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

ω = omega


"""
    mortality_vector(vec; start_age=0)

A convienience constructor to create an OffsetArray which has the array indexed by attained age rather than always starting from `1`. The package and JuliaActuary ecosystem assume that the rates are indexed by attained age and this allows transformation of tables without a direct dependency on **OffsetArrays.jl**.

Equivalent to doing:
```
using OffsetArrays
OffsetArray(vec,start_age-1)
```
"""

mortality_vector(vec; start_age=0) = return OffsetArray(vec,start_age-1)

"""
    reverse_lookup(table_dict,attained_age,rate;comparison)

Checks every table in the `key=>table` dictionary passed as `table_dict` for the given `rate` at `attained_age`. Will return a sorted array of tuples: 

- The key of the pair is the key of the table

- The first value of the tuple is the result of the comparison (for `isequals` will be `true` or `false`)
- The second value of the tuple is whether the table was a select or ultimate structure, if any
- The third value of the tuple is the issue_age, if table was select
"""
function reverse_lookup(tables,attained_age,rate;comparison=isequal)
    for (k1,v1) in tables
        reverse_lookup(table,attained_age,rate;comparison=comparison)
    end
end

function reverse_lookup(table::SelectUltimateTable{S,U},attained_age,rate;comparison=comparison) where {S,U}
    return (result=true,select=true,issue_age=4)
end

function reverse_lookup(table::UltimateTable{U},attained_age,rate;comparison=comparison) where {U}
    try
        return (result = comparison(rate,table.ultimate[attained_age]), select=false,issue_age=nothing)
    catch
        return nothing
    end
end

function reverse_lookup(vec::AbstractArray,attained_age,rate;comparison=comparison) where {U}
    try 
        result = comparison(rate,vec[attained_age])
    catch e
        result = nothing
    end
    return (result=result,select=nothing,issue_age=nothing)
end