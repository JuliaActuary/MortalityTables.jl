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

mort_max_dur = 121
mort_max_issue_age = 120
