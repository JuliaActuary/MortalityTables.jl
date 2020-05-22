"""
    get_SOA_table(id)

Given the id of a `mort.SOA.org` table, grab it and return it as a `MortalityTable`.

!! Remember that not all tables have been tested to work.
"""
function get_SOA_table(id::Int)
    path = "https://mort.soa.org/Export.aspx?Type=xml&TableIdentity=$id"
    r = HTTP.request("GET", path)

    # Why skip the first three bytes of the response?

    # From https://docs.python.org/3/library/codecs.html
    # To increase the reliability with which a UTF-8 encoding can be detected,
    # Microsoft invented a variant of UTF-8 (that Python 2.5 calls "utf-8-sig")
    # for its Notepad program: Before any of the Unicode characters is written
    #to the file, a UTF-8 encoded BOM (which looks like this as a byte sequence:
    # 0xef, 0xbb, 0xbf) is written.
    xml = getXML(String(r.body[4:end]))
    return XTbML_Table_To_MortalityTable(parseXTbMLTable(xml, path))

end

"""
    get_SOA_table!(dict,id)

Will lookup the given mortality table and add it to the given dict, with the name of
the table acting as the added key in the dictionary.

This modifies the given `dict` (as is indicated by the conventional `!` at the end
of the function name).

!! Remember that not all tables have been tested to work.
"""
function get_SOA_table!(dict,id::Int)
    path = "https://mort.soa.org/Export.aspx?Type=xml&TableIdentity=$id"
    r = HTTP.request("GET", path)

    # Why skip the first three bytes of the response?

    # From https://docs.python.org/3/library/codecs.html
    # To increase the reliability with which a UTF-8 encoding can be detected,
    # Microsoft invented a variant of UTF-8 (that Python 2.5 calls "utf-8-sig")
    # for its Notepad program: Before any of the Unicode characters is written
    #to the file, a UTF-8 encoded BOM (which looks like this as a byte sequence:
    # 0xef, 0xbb, 0xbf) is written.
    xml = getXML(String(r.body[4:end]))
    tbl = XTbML_Table_To_MortalityTable(parseXTbMLTable(xml, path))
    merge!(dict, Dict(tbl.d.name => tbl))

end
