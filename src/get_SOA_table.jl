"""
    get_SOA_table(id)
    get_SOA_table(table_name)

Given the id or name of a `mort.SOA.org` table, grab it and return it as a `MortalityTable`.

!!! Remember that not all tables have been tested to work.
"""
function get_SOA_table(id::Int)
    readXTbML(joinpath(artifact"mort.soa.org","t$id.xml"))
end

function get_SOA_table(table_name::String;source_map=table_source_map)
    entry = findfirst(x-> x.name == table_name, source_map)
    if entry === nothing
        search_method = StringDistances.Partial(StringDistances.Levenshtein())
        suggestion,_ = StringDistances.findnearest(table_name,[x.name for x in table_source_map],search_method)
        throw(ArgumentError("table name \"table_name\" not found in table set; " *
            "most similar available name is: \"$suggestion\""))
    end
    readXTbML(joinpath(artifact"mort.soa.org","t$(source_map[entry].id).xml"))
end

"""
    table(id)
    table(name)

Given the id or name of a `mort.SOA.org` table, grab it and return it as a `MortalityTable`.

!!! Remember that not all tables have been tested to work.
"""
table = get_SOA_table