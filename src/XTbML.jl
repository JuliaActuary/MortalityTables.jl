include("MortalityTable.jl")

function open_and_read(path)
    s = open(path) do file
        read(file, String)
    end
end

function getXML(open_file)

    return xml = XMLDict.xml_dict(open_file)

end

# get potentially missing value out of dict
function getVal(dict, key)
    try
        return val = parse(Float64, dict[key])
    catch y
        if isa(y, KeyError)
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

function parseXTbMLTable(x, path)
    md = x["XTbML"]["ContentClassification"]
    name = get(md, "TableName", nothing) |> strip
    id = get(md, "TableIdentity", nothing) |> strip
    provider = get(md, "ProviderName", nothing) |> strip
    reference = get(md, "TableReference", nothing) |> strip
    description = get(md, "TableDescription", nothing) |> strip
    comments = get(md, "Comments", nothing) |> strip
    source_path = path
    d = TableMetaData(
        name = name,
        id = id,
        provider = provider,
        reference = reference,
        description = description,
        comments = comments,
        source_path = source_path,
    )
    tbl = XTbMLTable(
        DataStructures.DefaultOrderedDict(missing),
        DataStructures.DefaultOrderedDict(missing),
        d,
    )

    if isa(x["XTbML"]["Table"], Vector)
        # for a select and ultimate table, will have multiple tables
        # parsed into a vector of tables
        for ai in x["XTbML"]["Table"][1]["Values"]["Axis"]
            issue_age = parse(Int, ai[:t])
            tbl.select[issue_age] = DataStructures.DefaultOrderedDict(missing)
            for aj in ai["Axis"]["Y"]
                duration = parse(Int, aj[:t])
                rate = getVal(aj, "")
                if !ismissing(rate)
                    tbl.select[issue_age][duration] = rate
                end
            end
        end

        for ai in x["XTbML"]["Table"][2]["Values"]["Axis"]["Y"]
            age = parse(Int, ai[:t])
            rate = getVal(ai, "")
            if !ismissing(rate)
                tbl.ultimate[age] = rate
            end
        end
    elseif isa(x["XTbML"]["Table"], DataStructures.OrderedDict)
        # a table without select period will just have one set of values, which
        # are loaded into the OrderedDict

        for ai in x["XTbML"]["Table"]["Values"]["Axis"]["Y"]
            age = parse(Int, ai[:t])
            rate = getVal(ai, "")
            if !ismissing(rate)
                tbl.ultimate[age] = rate
            end
        end
    else
        error("don't know how to handle table: " * name)
    end

    return tbl
end

"""
    q_select(table::XTbMLTable, issueAge, duration)
Given a mortality table, an issue age, and a duration, returns the appropriate select or ultimate qx.
"""
function q_select(table::XTbMLTable, issueAge::Int, duration::Int)
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
        q = table.ultimate[issueAge+duration-1]
    end
    return q
end

"""
    q_ulitmate(table::XTbMLTable, age)
Given a mortality table and an age returns the appropriate ultimate qx.
"""
function q_ultimate(table::XTbMLTable, age)
    return table.ultimate[age]
end

function XTbML_Table_To_MortalityTable(tbl::XTbMLTable)
    ult_α, ult_ω = extrema(keys(tbl.ultimate))
    ult = UltimateMortality([tbl.ultimate[age] for age = ult_α:ult_ω], ult_α)

    if length(tbl.select) > 0
        sel_α, sel_ω = extrema(keys(tbl.select))
        sel_end_dur = maximum([length(tbl.select[age]) for age = sel_α:sel_ω])

        select_array = Array{Any}(undef, length(sel_α:sel_ω), sel_end_dur)
        fill!(select_array, missing)

        for (iss_age, rate_vec) in tbl.select
            for (dur, rate) in rate_vec
                select_array[iss_age-sel_α+1, dur] = rate
            end
        end



        sel = SelectMortality(identity.(select_array), ult, sel_α)

        return MortalityTable(sel, ult, tbl.d)
    else
        return MortalityTable(ult, tbl.d)
    end
end

function readXTbML(path)
    path
    x = open_and_read(path) |> getXML
    XTbML_Table_To_MortalityTable(parseXTbMLTable(x, path))
end


# Load Available Tables ###

"""
    tables(dir=nothing)

    Loads the XtbML (the SOA XML data format for mortality tables) stored in the
    given path. If no path is specified, will load the packages in the
    MortalityTables package directory. To see where your system keeps packages,
    run `DEPOT_PATH` from a Julia REPL.
"""
function tables(dir = nothing)
    if isnothing(dir)
        table_dir = joinpath(dirname(pathof(MortalityTables)), "tables", "SOA")
    else
        table_dir = dir
    end
    tables = []
    @info "Loading built-in Mortality Tables..."
    for (root, dirs, files) in walkdir(table_dir)
        transducer = Filter(x -> basename(x)[end-3:end] == ".xml") |> Map(x ->  readXTbML(joinpath(root, x)) )
        tables = tcopy(transducer,files  )
    end
    # return tables
    return Dict(tbl.d.name => tbl for tbl in tables if ~isnothing(tbl))
end