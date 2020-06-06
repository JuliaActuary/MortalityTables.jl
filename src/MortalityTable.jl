include("MetaData.jl")


"""
Given an ultimate vector, will create a dictionary that is
indexed by issue age and will return `missing` `if the age is
not available.
"""
function UltimateMortality(v::Array{<:Real,1}, start_age = 0)
    return OffsetArray(v,start_age - 1)
end

"""
Given an 2D array, will create a an array that is indexed by issue age cotaining an array
which is then indexed by attained age.
"""
function SelectMortality(select, ultimate, start_age = 0)

    # iterate down the rows (issue ages)
    vs = map(enumerate(eachrow(select))) do (i, r)
        end_age = start_age + (i - 1) + (length(r) - 1)
        OffsetArray([r ; ultimate[end_age+1:end]],(start_age - 1) + (i - 1))
    end

    return OffsetArray(vs,start_age - 1)
end



"""
    MortalityTable

    A struct that holds a select (two-dimensional) and ultimate (vector) rates,
        along with MetaData associated with the table.
"""
abstract type MortalityTable end

struct SelectUltimateMortalityTable{S,U} <: MortalityTable
    select::S
    ultimate::U
    d::TableMetaData
end

struct UltimateMortalityTable{U} <: MortalityTable
    ultimate::U
    d::TableMetaData
end
function MortalityTable(select,ultimate,d::TableMetaData)
    return SelectUltimateMortalityTable(select, ultimate, d)
end

function MortalityTable(ultimate, d::TableMetaData)
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



function survivorship(v,to_age)
    return survivorship(v, firstindex(v), to_age)
end

function survivorship(v::T,from_age::Int,to_age::Int) where {T<:AbstractArray}
    if from_age == to_age
        return 1.0
    else
        return reduce(*,
            1 .- v[from_age:(to_age-1)],
            init=1.0
            )
    end
end

cumulative_decrement(v,to_age) = 1 .- survivorship(v,to_age) 
cumulative_decrement(v,from_age,to_age) = 1 .- survivorship(v,from_age,to_age) 

"""
    omega(x)
    ω(x)

Returns the last index of the given vector.
"""
function omega(x)
    return lastindex(x)
end

ω = omega