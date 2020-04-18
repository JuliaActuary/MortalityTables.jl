using DataStructures


include("MetaData.jl")


## refactor see https://github.com/JuliaActuary/MortalityTables.jl/issues/23
"""
MortalityVector
"""

struct MortalityVector
    q
end

"""
"""

abstract type MortalityDict end

struct SelectMortality <: MortalityDict
    v
end

struct UltimateMortality <: MortalityDict
    v
end

"""
Given an ultimate vector, will create a dictionary that is
indexed by issue age and will return `missing` `if the age is
not available.
"""
function UltimateMortality(v::Array{<:Real,1},start_age=0)
    d = DefaultDict(missing)
    for (i,q) in enumerate(v)
        d[start_age + i - 1] = MortalityVector(v[i:end])
    end
    return UltimateMortality(d)
end

"""
Given an 2D array, will create a dictionary that is
indexed by issue age and will return `missing` `if the age is
not available.
"""

function SelectMortality(select,ultimate::UltimateMortality,start_age=0)
    d = DefaultDict(missing)
    last_select_age = size(select,2) - 1  + size(select,1) - 1 + start_age

    # get the end of the table that would apply to the last select attained age
    last_ult_age = ω(ultimate,start_age + size(select,1)-1)

    # iterate down the rows (issue ages)
    for i in 1:size(select,1)
        iss_age = start_age + i - 1
        ult_age_start = length(select[i,1:end]) + iss_age
        select_qs = select[i,1:end]

        if ult_age_start >= maximum(keys(ultimate.v))
            d[iss_age] = select_qs |> MortalityVector

        else # use ultimate rates if available
            last_age = ω(ultimate,ult_age_start)
            last_dur = last_age - ult_age_start + 1
            start_dur = 1
            ult_qs = q(ultimate,ult_age_start,start_dur:last_dur)
            d[iss_age] = vcat(select_qs,ult_qs) |> MortalityVector
        end


    end

    return SelectMortality(d)
end



"""
    MortalityTable

    A struct that holds a select (two-dimensional) and ultimate (vector) rates,
        along with MetaData associated with the table.
"""
abstract type MortalityTable end

struct SelectUltimateMortalityTable <: MortalityTable
    select::SelectMortality
    ultimate::UltimateMortality
    d::TableMetaData
end

struct UltimateMortalityTable <: MortalityTable
    ultimate::UltimateMortality
    d::TableMetaData
end
function MortalityTable(select::SelectMortality,ultimate::UltimateMortality,d::TableMetaData)
    return SelectUltimateMortalityTable(select,ultimate,d)
end

function MortalityTable(ultimate::UltimateMortality,d::TableMetaData)
    # sel_α, sel_ω = extrema(keys(ultimate.v))
    # create a dummy select table which has the ultimate rate for the first
    # duration. From there, the normal SelectMortality constructor can take over

    # select = [q(ultimate,age,1) for age in sel_α:sel_ω ]
    # select = SelectMortality(select,ultimate,sel_α)
    return UltimateMortalityTable(ultimate,d)
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


##################################
### Basic Single Life Mortality ##
##################################

@doc raw"""
The probability that a life with given `issue_age` and currently in its nth
`duration`dies survives to at least `duration` + `time`. If given select
mortality, will be based on select rates.

Equivalant actuarial notation:
``$_tp_{(x)+s}$``, or the probability that a life aged `x + s` who was select
at age `x` survives to at least age `x+s+t`
"""
function p(table::MortalityDict,issue_age,duration,time::Int)
    prod(1.0 .- q(table,issue_age,duration:(duration+time-1)))
end
function p(table::MortalityDict,issue_age,duration,time)
    error("If you use non-integer times, you need to specify a \n
    distribution of deaths assumption (e.g. `Balducci()`, \n
    `Constant()`, or `Uniform()` as the last argument to your \n
    function call.")
end

@doc raw"""
the probability that a life aged `issue_age` + `duration` - 1
survives one additional timepoint

Equivalent actuarial notation:
``$p_x$`` , or
"""
function p(table::MortalityDict,issue_age,duration)
    return 1.0 - q(table,issue_age,duration)
end

function p(table::UltimateMortalityTable,issue_age,duration::Int)
    return p(table.ultimate,issue_age,duration)
end

function p(table::UltimateMortalityTable,issue_age,duration)
    error("If you use non-integer times, you need to specify a \n
    distribution of deaths assumption (e.g. `Balducci()`, \n
    `Constant()`, or `Uniform()` as the last argument to your \n
    function call.")
end

@doc raw"""
The probability that a life with given `issue_age` and currently in its nth
`duration`dies by at least `duration` + `time`. If given select mortality,
will be based on select rates.

Equivalent actuarial notation:
``$p_{(x)+s}$``  or the probability that a life aged `x + s` who was select
at age `x` dies by least age `x+s+t`
"""
function q(table::MortalityDict,issue_age,duration,time)
    1.0 - p(table::MortalityDict,issue_age,duration,time)
end

function q(table::UltimateMortalityTable,issue_age,duration,time)
    return q(table.ultimate,issue_age,duration,time)
end


function q(m::MortalityDict,issue_age::Int,duration)
    mv = m.v[issue_age]
    if ismissing(mv)
        return mv
    else
        return mv.q[duration]
    end
end

function q(m::UltimateMortality,issue_age::Int)
    mv = m.v[issue_age]
    if ismissing(mv)
        return mv
    else
        return mv.q[1]
    end
end

function q(m::UltimateMortalityTable,issue_age,duration)
    return q(m.ultimate,issue_age,duration)
end

function q(m::MortalityDict,issue_age::AbstractArray,duration)
    return [q(m,ia,duration) for ia in issue_age]
end

function q(m::UltimateMortalityTable,issue_age::AbstractArray,duration)
    return q(m.ultimate,issue_age,duration)
end

function omega(m::MortalityDict,issue_age)
    mv = m.v[issue_age]
    if ismissing(mv)
        return mv
    else
        return issue_age + length(m.v[issue_age].q) - 1
    end
end

function omega(m::UltimateMortalityTable,issue_age)
    return omega(m.ultimate,issue_age)
end

ω = omega
