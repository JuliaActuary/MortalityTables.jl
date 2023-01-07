"""
    get_SOA_table(id)
    get_SOA_table(table_name)

Given the id or name of a `mort.SOA.org` table, grab it and return it as a `MortalityTable`.

!!! Remember that not all tables have been tested to work.
"""
function get_SOA_table(id::Int)
    readXTbML(joinpath(artifact"mort.soa.org", "t$id.xml"))
end

function get_SOA_table2(id::Int)
    path = joinpath(artifact"mort.soa.org", "t$id.xml")
    leading_bytes = read(path,3)
    skipbytes = leading_bytes == [0xef, 0xbb, 0xbf]
            # Why skip the first three bytes of the response?

        # From https://docs.python.org/3/library/codecs.html
        # To increase the reliability with which a UTF-8 encoding can be detected,
        # Microsoft invented a variant of UTF-8 (that Python 2.5 calls "utf-8-sig")
        # for its Notepad program: Before any of the Unicode characters is written
        # to the file, a UTF-8 encoded BOM (which looks like this as a byte sequence:
        # 0xef, 0xbb, 0xbf) is written.

    x = open(path,"r") do f
        skipbytes && skip(f,3)
        XML.XMLTokenIterator(f) |> XML.Document
    end

    # t = parseXTbMLTable2(x,path)
    # XTbML_Table_To_MortalityTable(t)

end

function get_SOA_table(table_name::String; source_map = table_source_map)
    entry = get(source_map, table_name, nothing)
    if entry === nothing
        search_method = StringDistances.Partial(StringDistances.Levenshtein())
        suggestion, _ = StringDistances.findnearest(table_name,collect(keys(source_map)), search_method)
        throw(ArgumentError("table name \"table_name\" not found in table set; " *
                            "most similar available name is: \"$suggestion\""))
    end
    readXTbML(joinpath(artifact"mort.soa.org", "t$entry.xml"))
end

"""
    table(id)
    table(name)

Given the id or name of a `mort.SOA.org` table, grab it and return it as a `MortalityTable`.

!!! Remember that not all tables have been tested to work.
"""
table = get_SOA_table