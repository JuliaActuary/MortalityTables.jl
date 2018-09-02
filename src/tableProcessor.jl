import DataStructures
cso_vbt_2001 = DataStructures.OrderedDict(
    "set" => r"^(\d*\s\w{0,3})",
    "structure" => r"((Select|Ulitmate).*)(?=-)",
    "risk" => r"(.*)\s-" ,
    "sex" => r"(\w*)\s",
    "smoker" => r"(\w*)\,",
    "basis" => r"(\w*)")


#Given a Dict comprising the name and regex rule, will parse a string into the components
# e.g. "2001 CSO Preferred Select and Ultimate - Male Nonsmoker, ANB" and relevant ruleset
# becomes  ["2001 CSO", "Select and Ultimate", "Preferred", "Male", "Nonsmoker", "ANB"]
# and ["set", "structure", "risk", "sex", "smoker", "basis"]

function nameProcessor(name,rules)
    name = strip(name) # remove leading/trailing spaces
    ident_set = []
    parsed = []
    for (i,rule) in rules
        if occursin(rule,name)
            m = match(rule,name)
            name = strip(replace(name,rule => ""))

            append!(parsed,[strip(m[1])])
            append!(ident_set,[i])
        end
    end
    return ident_set, parsed
end
