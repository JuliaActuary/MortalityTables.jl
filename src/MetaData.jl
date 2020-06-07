"""
    TableMetaData(kwargs...)

Has the following fields, which default to `nothing` if not specified with a keyword:
- `name` - a name for the table
- `id` - if a mort.SOA.org sourced table, will be the identifying table ID
- `provider` - Where the rates came from
- `reference` - Source for more info on table
- `content_type`
- `description`
- `comments`
- `source_path`

When you call a `MortalityTable` interactively, it will nicely print this summary infomration.

# Example content from mort.SOA.org:

- **Table Identity**: 1076
- **Provider Domain**: actuary.org
- **Provider Name**: American Academy of Actuaries
- **Table Reference**: Tillinghast, “American Council of Life Insurers: Preferred Version of 2001 CSO Mortality Tables”, ...
- **Content Type**: CSO/CET
- **Table Name**: 2001 CSO Super Preferred Select and Ultimate - Male Nonsmoker, ANB
- **Table Description**: 2001 Commissioners Standard Ordinary (CSO) Super Preferred Select and Ultimate Table – Male Nonsmoker. Basis: Age Nearest Birthday. Minimum Select Age: 0. Maximum Select Age: 99. Minimum Ultimate Age: 16. Maximum Ultimate Age: 120
- **Comments**: Study Data: A preferred version of the 2001 Commissioners Standard Ordinary (CSO) table ...

And the `source_path` would be: [https://mort.soa.org/ViewTable.aspx?&TableIdentity=1076](https://mort.soa.org/ViewTable.aspx?&TableIdentity=1076)

# Example usage:
```
julia-repl> TableMetaData(name="My Table Name")
TableMetaData("My Table Name", nothing, nothing, nothing, nothing, nothing, nothing, nothing)
```

"""
Base.@kwdef struct TableMetaData    
    name::Union{Nothing,String} = nothing
    id::Union{Nothing,String} = nothing
    provider::Union{Nothing,String} = nothing
    reference::Union{Nothing,String} = nothing
    content_type::Union{Nothing,String} = nothing
    description::Union{Nothing,String} = nothing
    comments::Union{Nothing,String} = nothing
    source_path::Union{Nothing,String} = nothing
end