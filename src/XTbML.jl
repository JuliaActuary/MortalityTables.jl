

function open_and_read(path)
    bytes = read(path)
    if bytes[1:3] == [0xef, 0xbb, 0xbf]
        # Why skip the first three bytes of the response?

        # From https://docs.python.org/3/library/codecs.html
        # To increase the reliability with which a UTF-8 encoding can be detected,
        # Microsoft invented a variant of UTF-8 (that Python 2.5 calls "utf-8-sig")
        # for its Notepad program: Before any of the Unicode characters is written
        # to the file, a UTF-8 encoded BOM (which looks like this as a byte sequence:
        # 0xef, 0xbb, 0xbf) is written.
        return String(bytes[4:end])
    else
        return String(bytes)
    end
end

function getXML(open_file)

    return xml = XMLDict.xml_dict(open_file)

end

# get potentially missing value out of dict
function get_and_parse(dict, key)
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

struct XTbMLTable{S,U}
    select::S
    ultimate::U
    d::TableMetaData
end

function parseXTbMLTable(x, path)
    md = x["XTbML"]["ContentClassification"]
    name = get(md, "TableName", nothing) |> strip
    content_type = get(get(md, "ContentType", nothing), "", nothing) |> strip
    id = get(md, "TableIdentity", nothing) |> strip
    provider = get(md, "ProviderName", nothing) |> strip
    reference = get(md, "TableReference", nothing) |> strip
    description = get(md, "TableDescription", nothing) |> strip
    comments = get(md, "Comments", nothing) |> strip
    source_path = path
    d = TableMetaData(
        name=name,
        id=id,
        provider=provider,
        reference=reference,
        content_type=content_type,
        description=description,
        comments=comments,
        source_path=source_path,
    )

    if isa(x["XTbML"]["Table"], Vector)
        # for a select and ultimate table, will have multiple tables
        # parsed into a vector of tables
        sel = map(x["XTbML"]["Table"][1]["Values"]["Axis"]) do ai
            (issue_age = parse(Int, ai[:t]),
                rates = [(duration = parse(Int, aj[:t]), rate = get_and_parse(aj, "")) for aj in ai["Axis"]["Y"] if !ismissing(get_and_parse(aj, ""))])
        end

        ult = map(x["XTbML"]["Table"][2]["Values"]["Axis"]["Y"]) do ai 
            (age  = parse(Int, ai[:t]), rate = get_and_parse(ai, ""),)
        end

    else
        # a table without select period will just have one set of values

        ult = map(x["XTbML"]["Table"]["Values"]["Axis"]["Y"]) do ai
            (age = parse(Int, ai[:t]), 
                rate = get_and_parse(ai, ""))
        end

        sel = nothing

    end

    tbl = XTbMLTable(
        sel,
        ult,
        d
    )

    return tbl
end

# XTbML labels every rate with its age or duration. Each value goes to its label, so a table
# whose labels skip (rates at grouped ages, or an empty cell inside a select row) keeps its
# values at the right ages, with `missing` where it gives no rate. The element type widens only
# when there are such gaps. A repeated label is ambiguous and throws. `first` is the label of
# the first position (durations start at 1; ages at the smallest label).
function _by_label(labels, values, what, table; first = minimum(labels))
    allunique(labels) || throw(
        ArgumentError(
            "XTbML table $(something(table.name, table.source_path, "(unnamed)")): $what repeat a " *
                "label, so their rates cannot be placed: $(labels)"
        )
    )
    n = maximum(labels) - first + 1
    labels == first:(first + n - 1) && return first, collect(values)
    placed = Vector{Union{Missing, eltype(values)}}(missing, n)
    for (label, value) in zip(labels, values)
        placed[label - first + 1] = value
    end
    return first, placed
end

function XTbML_Table_To_MortalityTable(tbl::XTbMLTable)
    start_age, ult_rates = _by_label([v.age for v in tbl.ultimate], [v.rate for v in tbl.ultimate], "the ultimate ages", tbl.d)
    ult = UltimateMortality(ult_rates, start_age = start_age)

    if !isnothing(tbl.select)
        rows = map(tbl.select) do (issue_age, rates)
            # empty cells were dropped when parsing: durations without a rate are `missing`
            _, select_rates = _by_label(
                [r.duration for r in rates], [r.rate for r in rates],
                "the select durations for issue age $issue_age", tbl.d; first = 1
            )
            return _select_row(issue_age, select_rates, ult)
        end
        first_issue_age, sel = _by_label([r.issue_age for r in tbl.select], rows, "the select issue ages", tbl.d)
        sel = OffsetArray(sel, first_issue_age - 1)

        return MortalityTable(sel, ult, metadata=tbl.d)
    else
        return MortalityTable(ult, metadata=tbl.d)
    end
end

function _read_xtbml(path)
    x = open_and_read(path) |> getXML
    XTbML_Table_To_MortalityTable(parseXTbMLTable(x, path))
end

# Tables parsed from disk are cached by path so that repeated lookups of the
# same table return the same object without re-parsing the file.
const _TABLE_CACHE = Dict{String,Any}()
const _CACHE_LOCK = ReentrantLock()

"""
    readXTbML(path)

Loads the [XtbML](https://mort.soa.org/About.aspx) (the SOA XML data format for mortality tables) stored at the given path and returns a `MortalityTable`.

The result is cached by `path`, so calling this twice with the same path returns the identical object.
"""
function readXTbML(path)
    lock(_CACHE_LOCK) do
        get!(() -> _read_xtbml(path), _TABLE_CACHE, path)
    end
end


# Load Available Tables ###

"""
    read_tables(dir=nothing)

Loads the [XtbML](https://mort.soa.org/About.aspx) (the SOA XML data format for mortality tables) stored in the given path. If no path is specified, will load the packages in the MortalityTables package directory. To see where your system keeps packages, run `DEPOT_PATH` from a Julia REPL.
"""
function read_tables(dir=nothing)
    if isnothing(dir)
        table_dir = artifact"mort.soa.org"
    else
        table_dir = dir
    end
    tables = []
    @info "Loading built-in Mortality Tables..."
    for (root, dirs, files) in walkdir(table_dir)
        for file in files
            if endswith(file,".xml") && !startswith(file,".")
                tbl =  readXTbML(joinpath(root,file))
                push!(tables,tbl)
            end
        end
    end
    return Dict(tbl.metadata.name => tbl for tbl in tables if ~isnothing(tbl))
end


# this is used to generate the table mapping in table_source_map.jl
function _write_available_tables()
        table_dir = artifact"mort.soa.org"
    tables = []
    @info "Loading built-in Mortality Tables..."
    for (root, dirs, files) in walkdir(table_dir)
        for file in files
            if endswith(file,".xml") && !startswith(file,".")
            x = open_and_read(joinpath(root,file)) |> XMLDict.xml_dict
            md = x["XTbML"]["ContentClassification"]
            name = get(md, "TableName", nothing) |> strip
            content_type = get(get(md, "ContentType", nothing), "", nothing) |> strip
            id = get(md, "TableIdentity", nothing) |> strip
            push!(tables,(source="mort.soa.org",name=name,id=parse(Int,id)))
            end
        end
    end
    return sort!(tables,by=last)
end