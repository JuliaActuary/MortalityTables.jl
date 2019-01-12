import DataStructures

include("XTbMl.jl")

if Int === Int64
    struct MortalityTable
        select::Array{Union{Missing, Float64},2}
        ultimate::Array{Union{Missing, Float64},1}
    end
else
    struct MortalityTable
        select::Array{Union{Missing, Float32},2}
        ultimate::Array{Union{Missing, Float32},1}
    end
end

"""
    qx(table::XTbMLTable, issueAge, duration)

Given a mortality table, an issue age, and a duration, returns the appropriate select or ultimate qx.
"""
function qx(table::XTbMLTable, issueAge::Int, duration::Int)

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

# Load Available Tables
table_dir = joinpath(dirname(pathof(MortalityTables)), "tables", "SOA")

function Tables()
    tables = Dict()
    for (root, dirs, files) in walkdir(table_dir)
        for file in files

            if basename(file)[end-3:end] == ".xml"
                tbl, name = loadXTbMLTable(joinpath(root,file))
                tables[strip(name)] = XtbML_Table_To_Matrix(tbl) #strip removes leading/trailing whitespace from the name
            end
        end
    end
    return tables
end



function XtbML_Table_To_Matrix(tbl::XTbMLTable)
    return MortalityTable(ones(2,2),[qx(tbl,age) for age=0:120])
    # return MortalityTable([qx(tbl,issue_age,dur) for issue_age=0:120,dur=1:121],[qx(tbl,age) for age=0:120])
end

### TABLE STRUCUTRE PARSING ###
# extract table varies-by characteristics using a regex rule

cso_vbt_2001 = DataStructures.OrderedDict(
    "set" => r"^(\d*\s\w{0,3})",
    "structure" => r"((Select|Ulitmate).*)(?=-)",
    "risk" => r"(.*)\s-" ,
    "sex" => r"(\w*)\s",
    "smoker" => r"(\w*)\,",
    "basis" => r"(\w*)")


#Given a Dict comprising the name and regex rule, will parse a string into the components
# e.g. "2001 CSO Preferred Select and Ultimate - Male Nonsmoker, ANB" and relevant ruleset
# becomes  ["2001 CSO", "Select and Ultimate", "Preferred", "Male", "Nonsmoker", "ANB"]
# and ["set", "structure", "risk", "sex", "smoker", "basis"]

function nameProcessor(name,rules)
    name = strip(name) # remove leading/trailing spaces
    ident_set = []
    parsed = []
    for (i,rule) in rules
        if occursin(rule,name)
            m = match(rule,name)
            name = strip(replace(name,rule => ""))

            append!(parsed,[strip(m[1])])
            append!(ident_set,[i])
        end
    end
    return ident_set, parsed
end
