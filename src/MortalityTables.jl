module MortalityTables
import XMLDict
import DataStructures

include("tableProcessor.jl")

"""
    qx(table::XTbMLTable, issueAge, duration)

Given a mortality table, an issue age, and a duration, returns the appropriate select or ultimate qx.
"""
function qx(table::XTbMLTable, issueAge, duration)

    if length(table.select) > 0
        q = table.select[issueAge][duration]
    else
        q = missing
    end

    if ismissing(q)
        q = table.ultimate[issueAge + duration - 1]
    end
    return q
end


"""
    qx(table::XTbMLTable, age)

Given a mortality table and an age returns the appropriate ultimate qx.
"""
function qx(table::XTbMLTable, age)
    return table.ultimate[age]
end

end # module

export Tables,loadXTbMLTable, qx
