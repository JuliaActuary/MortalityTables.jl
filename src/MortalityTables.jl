module MortalityTables
import XMLDict
import DataStructures

include("tableProcessor.jl")

"""
    qx(table::MortalityTable, issueAge::Int, duration::Int)

Given a mortality table, an issue age, and a duration, returns the appropriate select or ultimate qx.

If the table is a select and ultimate table, it will return the select mortality during the select period.
If the table does not have select rates, will just return the rate for the age at the given duration.
"""

function qx(table::MortalityTable, issueAge::Int, duration::Int)
    try
        table.select[issueAge+1,duration] # +1 because julia index starts with 1
    catch err
        if isa(err,BoundsError)
            missing
        else
            table.select[issueAge+1,duration]
        end
    end
end


"""
    qx(table::XTbMLTable, age)

Given a mortality table and an age returns the qx.

If the table is a select and ultimate table, it will return the select mortality during the ulitmate period because duration was not given.

If the table does not have select rates, will just return the rate for the given age.

"""
function qx(table::MortalityTable, age)
    return table.ultimate[age+1] # +1 because julia index starts with 1
end

end # module

export Tables,loadXTbMLTable, qx
