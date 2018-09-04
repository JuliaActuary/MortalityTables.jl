import DataStructures
import XMLDict

function getXML(path)
    s = open(path) do file
        read(file,String)
    end
    return xml = XMLDict.xml_dict(s)

end
# get potentially missing value out of dict
function getVal(dict,key)
    try
        return val = parse(Float64,dict[key])
    catch y
        if isa(y,KeyError)
            return val = missing
        else
            throw(y)
        end
    end
end

struct XTbMLTable
    select::DataStructures.DefaultOrderedDict
    ultimate::DataStructures.DefaultOrderedDict
end

function loadXTbMLTable(path)

    x = getXML(path)
    name = x["XTbML"]["ContentClassification"]["TableName"]
    tbl = XTbMLTable(DataStructures.DefaultOrderedDict(missing),DataStructures.DefaultOrderedDict(missing))

    if isa(x["XTbML"]["Table"],Vector)
        # for a select and ultimate table, will have multiple tables
        # parsed into a vector of tables
        for ai in x["XTbML"]["Table"][1]["Values"]["Axis"]
            issue_age = parse(Int,ai[:t])
            tbl.select[issue_age] = DataStructures.DefaultOrderedDict(missing)
            for aj in ai["Axis"]["Y"]
                duration = parse(Int,aj[:t])
                rate = getVal(aj,"")
                if !ismissing(rate)
                    tbl.select[issue_age][duration]= rate
                end
            end
        end

        for ai in x["XTbML"]["Table"][2]["Values"]["Axis"]["Y"]
            age = parse(Int,ai[:t])
            rate = getVal(ai,"")
            if !ismissing(rate)
                tbl.ultimate[age]= rate
            end
        end
    elseif isa(x["XTbML"]["Table"],DataStructures.OrderedDict)
        # a table without select period will just have one set of values, which
        # are loaded into the OrderedDict

        for ai in x["XTbML"]["Table"]["Values"]["Axis"]["Y"]
            age = parse(Int,ai[:t])
            rate = getVal(ai,"")
            if !ismissing(rate)
                tbl.ultimate[age]= rate
            end
        end
    else
        error("don't know how to handle table: "  *  name)
    end

    return (tbl, name)
end

# tb = XTbMLTable(DataStructures.DefaultOrderedDict(missing),DataStructures.DefaultOrderedDict(missing))
# vbt = loadXTbMLTable(raw"C:\Users\alecl\AppData\Local\Julia-1.0.0\MortalityTables\src\tables\SOA\t17.xml")
# cso80 = getXML(raw"C:\Users\alecl\AppData\Local\Julia-1.0.0\MortalityTables\src\tables\SOA\t17.xml")
# cso01 = getXML(raw"C:\Users\alecl\AppData\Local\Julia-1.0.0\MortalityTables\src\tables\SOA\t1076.xml")
# cso80["XTbML"]["Table"]["Values"]["Axis"]["Y"]
# cso01["XTbML"]["Table"]
#
# typeof(cso80["XTbML"]["Table"])
# cso80
