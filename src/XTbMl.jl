import DataStructures
import XMLDict

include("MortalityTable.jl")

function open_and_read(path)
    s = open(path) do file
        read(file,String)
    end
end

function getXML(open_file)

    return xml = XMLDict.xml_dict(open_file)

end

# get potentially missing value out of dict
function getVal(dict,key)
    try
        return val = parse(Float64,dict[key])
    catch y
        if isa(y,KeyError)
            return val = missing
        else
            throw(y)
        end
    end
end

struct XTbMLTable
    select::DataStructures.DefaultOrderedDict
    ultimate::DataStructures.DefaultOrderedDict
    d::TableMetaData
end

function parseXTbMLTable(x,path)
    md = x["XTbML"]["ContentClassification"]
    name = get(md,"TableName",nothing)
    id = get(md,"TableIdentity",nothing)
    provider = get(md,"ProviderName",nothing)
    reference = get(md,"TableReference",nothing)
    description = get(md,"TableDescription",nothing)
    comments = get(md,"Comments",nothing)
    source_path = path
    d = TableMetaData(name=name,id=id,provider=provider,reference=reference,
        description=description,comments=comments,source_path=source_path)
    tbl = XTbMLTable(DataStructures.DefaultOrderedDict(missing),DataStructures.DefaultOrderedDict(missing),d)

    if isa(x["XTbML"]["Table"],Vector)
        # for a select and ultimate table, will have multiple tables
        # parsed into a vector of tables
        for ai in x["XTbML"]["Table"][1]["Values"]["Axis"]
            issue_age = parse(Int,ai[:t])
            tbl.select[issue_age] = DataStructures.DefaultOrderedDict(missing)
            for aj in ai["Axis"]["Y"]
                duration = parse(Int,aj[:t])
                rate = getVal(aj,"")
                if !ismissing(rate)
                    tbl.select[issue_age][duration]= rate
                end
            end
        end

        for ai in x["XTbML"]["Table"][2]["Values"]["Axis"]["Y"]
            age = parse(Int,ai[:t])
            rate = getVal(ai,"")
            if !ismissing(rate)
                tbl.ultimate[age]= rate
            end
        end
    elseif isa(x["XTbML"]["Table"],DataStructures.OrderedDict)
        # a table without select period will just have one set of values, which
        # are loaded into the OrderedDict

        for ai in x["XTbML"]["Table"]["Values"]["Axis"]["Y"]
            age = parse(Int,ai[:t])
            rate = getVal(ai,"")
            if !ismissing(rate)
                tbl.ultimate[age]= rate
            end
        end
    else
        error("don't know how to handle table: "  *  name)
    end

    return tbl
end

"""
    qx(table::XTbMLTable, issueAge, duration)
Given a mortality table, an issue age, and a duration, returns the appropriate select or ultimate qx.
"""
function qx(table::XTbMLTable, issueAge::Int, duration::Int)
    if length(table.select) > 0
        s = table.select[issueAge]
        if ismissing(s)
            q = missing
        else
            q = s[duration]
        end
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

function XTbML_Table_To_MortalityTable(tbl::XTbMLTable)

    return MortalityTable(
    OffsetArray([qx(tbl,issue_age,dur) for issue_age=0:mort_max_issue_age,dur=1:mort_max_dur],-1,0),
    OffsetArray([qx(tbl,age) for age=0:mort_max_issue_age],-1),
    tbl.d)
end

function readXTbML(path)
    x = open_and_read(path) |> getXML
    XTbML_Table_To_MortalityTable(parseXTbMLTable(x,path))
end


# Load Available Tables ###

"""
    tables(dir=nothing)

    Loads the XtbML (the SOA XML data format for mortality tables) stored in the
    given path. If no path is specified, will load the packages in the
    MortalityTables package directory. To see where your system keeps packages,
    run `DEPOT_PATH` from a Julia REPL.
"""
function tables(dir=nothing)
    if isnothing(dir)
        table_dir = joinpath(dirname(pathof(MortalityTables)), "tables", "SOA")
    else
        table_dir = dir
    end

    tables = Dict()
    for (root, dirs, files) in walkdir(table_dir)
        for file in files

            if basename(file)[end-3:end] == ".xml"
                tbl = readXTbML(joinpath(root,file))
                tables[tbl.d.name] = tbl
            end
        end
    end
    return tables
end
