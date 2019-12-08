using OffsetArrays


include("MetaData.jl")

"""
`MortalityMatrix` is a 2-dimensional array with duration across the columns
and issue age down the rows. Issue age must begin at age 0 and duration must
begin at 1. The index of the table is [`issue_age`,`duration`].

For a select table, duration will run across the row until the end of the rates (this is different than
how tables are typcially displayed given limited width constraints).
"""
const MortalityMatrix= Union{OffsetArray{Union{Missing, Float64},2,Array{Union{Missing, Float64},2}},
                                OffsetArray{Float64,2,Array{Float64,2}}}

const MortalityArray=  Union{OffsetArray{Union{Missing, Float64},1,Array{Union{Missing, Float64},1}},
                                OffsetArrays.OffsetArray{Float64,1,Array{Float64,1}}}

"""
    ultimate_vector_to_matrix(v,start_age=0)

Given a flat vector, turn it into a matrix to fit a select format.
This enables the same call for select mortality to work on an ultimate-only
    table. """
function ultimate_vector_to_matrix(v,start_age=0)
    age_offset = start_age - 1
    v = OffsetArray(v,start_age - 1)
    function inbounds(v,i)
        indices = axes(v)[1]
         first(indices) <= i <= last(indices)
    end
    OffsetArray([inbounds(v,iss_age + dur - 1) ? v[iss_age + dur - 1] : missing for dur=1:mort_max_dur, iss_age=0:mort_max_issue_age],-1,0)
end

function normalized_vector(v,start_age=0)
    age_offset = start_age - 1
    v = OffsetArray(v,start_age - 1)
    function inbounds(v,i)
        indices = axes(v)[1]
        first(indices) <= i <= last(indices)
    end

    OffsetArray([inbounds(v,age) ? v[age] : missing for age=0:mort_max_issue_age],-1)
end
"""
    MortalityTable

    A struct that holds a select (two-dimensional) and ultimate (vector) rates,
        along with MetaData associated with the table.
"""
struct MortalityTable
    select::MortalityMatrix
    ultimate::MortalityArray
    d::TableMetaData
end

"""
    MortalityTable(vector,start_age,::MetaData)

    A constructor that will convert a 1-d array of rates into a MortalityTable
    object. `start_age` is the first age for which the rates apply.
"""
function MortalityTable(ultimate::Array{Float64,1},start_age,name=TableMetaData())
    MortalityTable(
    ultimate_vector_to_matrix(ultimate,start_age),
    normalized_vector(ultimate,start_age),
    name)
end


Base.show(io::IO, ::MIME"text/plain", mt::MortalityTable) =
    print(io,
    """
    MortalityTable:
       Name:
           $(mt.d.name)
       Provider:
           $(mt.d.provider)
       mort.SOA.org ID:
           $(mt.d.id)
       mort.SOA.org link:
           https://mort.soa.org/ViewTable.aspx?&TableIdentity=$(mt.d.id)
       Description:
           $(mt.d.description)
    """)

"""
Return a vector/array of mortality rates. Two arguments returns select rates,
while one argument returns ultimate rates.
"""
function Base.getindex(mt::MortalityTable,ages,durs)
    return mt.select[ages,durs]
end

function Base.getindex(mt::MortalityTable,ages)
    return mt.ultimate[ages]
end


##################################
### Basic Single Life Mortality ##
##################################

"""
ₜp₍ₓ₎₊ₛ , or the probability that a life aged `x + s` who was select
at age `x` survives to at least age `x+s+t`
"""
function p(table::MortalityTable,x,s,t)
    prod(1.0 .- table.select[x,s+1:s+t])
end

"""
ₜpₓ , or the probability that a life aged `x` survives to at least age `t`
"""
function p(table::MortalityTable,x,t)
    prod(1.0 .- table.ultimate[x:x+t-1])
end

"""
pₓ , or the probability that a life aged `x` survives through age `x+1`
"""
function p(table::MortalityTable,x)
    1.0 - q(table,x)
end

"""
ₜq₍ₓ₎₊ₛ , or the probability that a life aged `x + s` who was select
at age `x` dies by least age `x+s+t`
"""
function q(table::MortalityTable,x,s,t)
    1.0 - p(table,x,s,t)
end

"""
ₜqₓ , or the probability that a life aged `x` dies by age `x+t`
"""
function q(table::MortalityTable,x,t)
    1.0 - p(table,x,t)
end

"""
qₓ , or the probability that a life aged `x` dies by age `x+1`
"""
function q(table::MortalityTable,x)
    table.ultimate[x]
end

"""
`qx` is a convenience function that allows you to get the rate at a given `age`.
If wanting select/ultimate rates, specifiy the `duration` and `age` should be the issue age.
"""
function qx(table::MortalityTable,age)
    q(table,age)
end

function qx(table::MortalityTable,age,duration)
    table.select[age,duration]
end
