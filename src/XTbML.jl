

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

### XTbML node helpers

_by_tag(node, name) = filter(e -> XML.tag(e) == name, XML.elements(node))
_child(node, name) = only(_by_tag(node, name))

# Text of an element's Text/CData children. Attributes are allowed, so
# `<ContentType tc="85">CSO/CET</ContentType>` works; XML.jl's `simple_value`
# rejects elements with attributes. The `replace` is XML 1.0 §2.11 line-end
# normalization, which libxml2 used to do for us and XML.jl 0.4.6 does not yet
# (merged upstream in JuliaData/XML.jl #133 and #136, not yet released).
_content(e) = replace(
    join(XML.value(c) for c in XML.children(e) if XML.nodetype(c) in (XML.Text, XML.CData)),
    "\r\n" => "\n", "\r" => "\n")

# Metadata field: `nothing` when absent, "" when present but empty.
_text(parent, name) = (es = _by_tag(parent, name); isempty(es) ? nothing : _content(first(es)))

# A rate cell: `<Y t="3">0.00260</Y>` is a rate; an empty or blank `<Y t="3"></Y>` is `missing`.
function _rate(y)
    s = strip(_content(y))
    return isempty(s) ? missing : parse(Float64, s)
end

function _content_classification(root, path)
    md = _child(root, "ContentClassification")
    field(name) = (s = _text(md, name); s === nothing ? nothing : String(strip(s)))
    return TableMetaData(
        name = field("TableName"),
        id = field("TableIdentity"),
        provider = field("ProviderName"),
        reference = field("TableReference"),
        content_type = field("ContentType"),
        description = field("TableDescription"),
        comments = field("Comments"),
        source_path = path,
    )
end

"""
    parseXTbMLTable(str, path)

Parse the XTbML document in `str` (read from `path`, which is recorded in the metadata) into a
named tuple `(select, ultimate, metadata)`. `ultimate` is a vector of `(age, rate)`; `select` is
`nothing` for an ultimate-only table, or a vector of `(issue_age, rates)` where `rates` is a vector
of `(duration, rate)` for the defined durations.
"""
function parseXTbMLTable(str::AbstractString, path)
    root = _child(XML.parse(str, XML.Node), "XTbML")
    tables = _by_tag(root, "Table")
    ys(tbl) = _by_tag(_child(_child(tbl, "Values"), "Axis"), "Y")
    ult = [(age = parse(Int, y["t"]), rate = _rate(y)) for y in ys(tables[end])]
    # a select and ultimate table has two <Table>s: the select rates (by issue
    # age, then duration) followed by the ultimate rates
    sel = length(tables) == 1 ? nothing : map(_by_tag(_child(tables[1], "Values"), "Axis")) do ai
        rates = [(duration = parse(Int, y["t"]), rate = _rate(y)) for y in _by_tag(_child(ai, "Axis"), "Y")]
        (issue_age = parse(Int, ai["t"]), rates = filter(r -> !ismissing(r.rate), rates))
    end
    return (select = sel, ultimate = ult, metadata = _content_classification(root, path))
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

function XTbML_Table_To_MortalityTable(tbl)
    start_age, ult_rates = _by_label([v.age for v in tbl.ultimate], [v.rate for v in tbl.ultimate], "the ultimate ages", tbl.metadata)
    ult = UltimateMortality(ult_rates, start_age = start_age)

    if !isnothing(tbl.select)
        rows = map(tbl.select) do (issue_age, rates)
            # empty cells were dropped when parsing: durations without a rate are `missing`
            _, select_rates = _by_label(
                [r.duration for r in rates], [r.rate for r in rates],
                "the select durations for issue age $issue_age", tbl.metadata; first = 1
            )
            return _select_row(issue_age, select_rates, ult)
        end
        first_issue_age, sel = _by_label([r.issue_age for r in tbl.select], rows, "the select issue ages", tbl.metadata)
        sel = OffsetArray(sel, first_issue_age - 1)

        return MortalityTable(sel, ult, metadata=tbl.metadata)
    else
        return MortalityTable(ult, metadata=tbl.metadata)
    end
end

_read_xtbml(path) = XTbML_Table_To_MortalityTable(parseXTbMLTable(open_and_read(path), path))

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
                path = joinpath(root, file)
                doc = XML.parse(open_and_read(path), XML.Node)
                md = _content_classification(_child(doc, "XTbML"), path)
                push!(tables, (source="mort.soa.org", name=md.name, id=parse(Int, md.id)))
            end
        end
    end
    return sort!(tables,by=last)
end
