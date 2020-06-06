include("MetaData.jl")


"""
"""

abstract type MortalityRates end

struct SelectMortality{T} <: MortalityRates
    v::T
end

struct UltimateMortality{T} <: MortalityRates
    v::T
end

"""
Given an ultimate vector, will create a dictionary that is
indexed by issue age and will return `missing` `if the age is
not available.
"""
function UltimateMortality(v::Array{<:Real,1}, start_age = 0)
    return OffsetArray(v,start_age) |> UltimateMortality
end

"""
Given an 2D array, will create a an array that is indexed by issue age cotaining an array
which is then indexed by attained age.
"""
function SelectMortality(select, ultimate::UltimateMortality, start_age = 0)

    last_select_age = size(select, 2) - 1 + size(select, 1) - 1 + start_age

    # get the end of the table that would apply to the last select attained age
    last_ult_age = ω(ultimate, start_age + size(select, 1) - 1)

    # iterate down the rows (issue ages)
    vs = map(enumerate(eachrow(select))) do (i, r)
        end_age = start_age + i - 1 + length(r)
        OffsetArray([r ; ultimate[end_age+1:end]],start_age + i)
    end

    return OffsetArray(vs,start_age) |> SelectMortality
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
function MortalityTable(
    select::SelectMortality,
    ultimate::UltimateMortality,
    d::TableMetaData,
)
    return SelectUltimateMortalityTable(select, ultimate, d)
end

function MortalityTable(ultimate::UltimateMortality, d::TableMetaData)
    # sel_α, sel_ω = extrema(keys(ultimate.v))
    # create a dummy select table which has the ultimate rate for the first
    # duration. From there, the normal SelectMortality constructor can take over

    # select = [q(ultimate,age,1) for age in sel_α:sel_ω ]
    # select = SelectMortality(select,ultimate,sel_α)
    return UltimateMortalityTable(ultimate, d)
end


Base.show(io::IO, ::MIME"text/plain", mt::MortalityTable) = print(
    io,
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
    """,
)


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
function p(table::MortalityRates, issue_age, duration, time::Int)
    if time == 0 
        return 1.0
    else
        reduce(*,1.0 .- q(table, issue_age, duration:(duration+time-1)))
    end
end

function p(table::MortalityRates, issue_age, duration, time)
    throw(ArgumentError("time: $time - If you use non-integer time, you need to specify a \n
          distribution of deaths assumption (e.g. `Balducci()`, \n
          `Constant()`, or `Uniform()` as the last argument to your \n
          function call."))
end

@doc raw"""
the probability that a life aged `issue_age` + `duration` - 1
survives one additional timepoint

Equivalent actuarial notation:
``$p_x$`` , or
"""
function p(table::MortalityRates, issue_age, duration)
    return 1.0 .- q(table, issue_age, duration)
end

function p(table::UltimateMortalityTable, issue_age, duration::Int)
    return p(table.ultimate, issue_age, duration)
end

function p(table::UltimateMortalityTable, issue_age, duration)
    throw(ArgumentError("time: $time - If you use non-integer time, you need to specify a distribution of deaths assumption (e.g. `Balducci()`, Constant()`, or `Uniform()` as the last argument to your function call."))
end

@doc raw"""
The probability that a life with given `issue_age` and currently in its nth
`duration`dies by at least `duration` + `time`. If given select mortality,
will be based on select rates.

Equivalent actuarial notation:
``$p_{(x)+s}$``  or the probability that a life aged `x + s` who was select
at age `x` dies by least age `x+s+t`
"""
function q(table::MortalityRates, issue_age, duration, time::Int)
    1.0 - p(table::MortalityRates, issue_age, duration, time)
end
function q(table::MortalityRates, issue_age, duration, time)
    throw(ArgumentError("time: $time - If you use non-integer time, you need to specify a distribution of deaths assumption (e.g. `Balducci()`, `Constant()`, or `Uniform()` as the last argument to your \n
          function call."))

end

function q(table::UltimateMortalityTable, issue_age, duration, time::Int)
    return q(table.ultimate, issue_age, duration, time)
end

function q(table::UltimateMortalityTable, issue_age, duration, time)
    throw(ArgumentError("time: $time - If you use non-integer time, you need to specify a distribution of deaths assumption (e.g. `Balducci()`, `Constant()`, or `Uniform()` as the last argument to your 
          function call."))
end


function q(m::MortalityRates, issue_age::Int, duration)
    mv = m.v[issue_age + duration - 1]
    if ismissing(mv)
        return mv
    else
        return mv[duration]
    end
end

function q(m::UltimateMortality, issue_age::Int)
    mv = m.v[issue_age]
    if ismissing(mv)
        return mv
    else
        return mv[1]
    end
end

function q(m::UltimateMortalityTable, issue_age, duration)
    return q(m.ultimate, issue_age, duration)
end

function q(m::MortalityRates, issue_age::AbstractArray, duration)
    return [q(m, ia, dur) for ia in issue_age, dur in duration]
end

function q(m::UltimateMortalityTable, issue_age, duration)
    return q(m.ultimate, issue_age, duration)
end

function omega(m::SelectMortality, issue_age)
    return lastindex(m.v[issue_age])
end

function omega(m::UltimateMortalityTable)
    return omega(m.ultimate.v)
end

ω = omega
