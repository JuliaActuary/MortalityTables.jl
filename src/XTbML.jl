

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
    field(name) = _metadata_text(_text(md, name))
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

Parse the XTbML document in `str` (read from `path`, which is recorded in the metadata) into the
labeled rates `(select, ultimate, metadata)` that `_table_from_labels` assembles. `ultimate` is a
vector of `(age, rate)`; `select` is `nothing` for an ultimate-only table, or a vector of
`(issue_age, rates)` where `rates` is a vector of `(duration, rate)`, `missing` for an empty cell.
"""
function parseXTbMLTable(str::AbstractString, path)
    root = _child(XML.parse(str, XML.Node), "XTbML")
    metadata = _content_classification(root, path)
    # the select rates are by issue age, then duration
    sel_table, ult_table = _select_and_ultimate(_by_tag(root, "Table"), "<Table> elements", _xtbml_source(metadata))
    ys(tbl) = _by_tag(_child(_child(tbl, "Values"), "Axis"), "Y")
    ult = [(age = parse(Int, y["t"]), rate = _rate(y)) for y in ys(ult_table)]
    sel = isnothing(sel_table) ? nothing : map(_by_tag(_child(sel_table, "Values"), "Axis")) do ai
        # the durations are read before the issue age, so a select table without a duration axis
        # fails on that missing <Axis>, the error test/data/unsupported_tables.txt records
        rates = [(duration = parse(Int, y["t"]), rate = _rate(y)) for y in _by_tag(_child(ai, "Axis"), "Y")]
        (issue_age = parse(Int, ai["t"]), rates = rates)
    end
    return (select = sel, ultimate = ult, metadata = metadata)
end

# How errors name an XTbML table: by its name, or else by the file it was read from.
_xtbml_source(metadata) = "XTbML table $(something(metadata.name, metadata.source_path, "(unnamed)"))"

function _read_xtbml(path)
    tbl = parseXTbMLTable(open_and_read(path), path)
    return _table_from_labels(tbl, _xtbml_source(tbl.metadata))
end

# Tables parsed from disk are cached by path so that repeated lookups of the
# same table return the same object without re-parsing the file.
const _TABLE_CACHE = Dict{String,MortalityTable}()
const _CACHE_LOCK = ReentrantLock()

"""
    readXTbML(path)

Loads the [XTbML](https://mort.soa.org/About.aspx) (the SOA XML data format for mortality tables) stored at the given path and returns a `MortalityTable`.

The result is cached by `path`, so calling this twice with the same path returns the identical object.
"""
function readXTbML(path)
    lock(_CACHE_LOCK) do
        get!(() -> _read_xtbml(path), _TABLE_CACHE, path)
    end
end


# Load Available Tables ###

# The XTbML files anywhere under `dir`, skipping hidden files (such as macOS `._` files).
_xtbml_paths(dir) =
    [joinpath(root, file) for (root, _, files) in walkdir(dir) for file in files if endswith(file, ".xml") && !startswith(file, ".")]

"""
    read_tables(dir=nothing)

Loads the [XTbML](https://mort.soa.org/About.aspx) files (the SOA XML data format for mortality tables) stored in the given directory and returns a `Dict` of the tables by name. If no directory is given, it loads the tables bundled with the package.
"""
function read_tables(dir=nothing)
    table_dir = isnothing(dir) ? artifact"mort.soa.org" : dir
    @info "Loading built-in Mortality Tables..."
    tables = MortalityTable[readXTbML(path) for path in _xtbml_paths(table_dir)]
    return Dict(tbl.metadata.name => tbl for tbl in tables)
end


# this is used to generate the table mapping in table_source_map.jl
function _write_available_tables()
    @info "Loading built-in Mortality Tables..."
    tables = map(_xtbml_paths(artifact"mort.soa.org")) do path
        md = _content_classification(_child(XML.parse(open_and_read(path), XML.Node), "XTbML"), path)
        (source = "mort.soa.org", name = md.name, id = parse(Int, md.id))
    end
    return sort!(tables, by = last)
end
