module MortalityTablesCSVExt

using MortalityTables, CSV, OffsetArrays
import MortalityTables: MortalityTable, TableMetaData, UltimateMortality

""" 
    MortalityTable(CSV.File)

Read and parse a CSV file in the SOA's usual CSV format. You must import and use CSV.jl before calling this function.

# Examples

```julia-repl
julia> path = "path/to/table.csv"
julia> file = CSV.File(path,header=false) # no real header in the file
julia> MortalityTable(file)
MortalityTable (Insured Lives Mortality):
   Name:
       2015 VBT Female Non-Smoker RR50 ALB
   Fields: 
       (:select, :ultimate, :metadata)
   Provider:
       American Academy of Actuaries along with the Society of Actuaries
   mort.SOA.org ID:
       3209
   mort.SOA.org link:
       https://mort.soa.org/ViewTable.aspx?&TableIdentity=3209
   Description:
       2015 VBT Relative Risk Table - Female, 50% Non-Smoker, Age Last Birthday, Select
```

And in one line: 
```julia-repl
julia> MortalityTables.MortalityTable(CSV.File( "path/to/table.csv",header=false))
MortalityTable (Insured Lives Mortality):
   Name:
       2015 VBT Female Non-Smoker RR50 ALB
   Fields: 
       (:select, :ultimate, :metadata)
   Provider:
       American Academy of Actuaries along with the Society of Actuaries
   mort.SOA.org ID:
       3209
   mort.SOA.org link:
       https://mort.soa.org/ViewTable.aspx?&TableIdentity=3209
   Description:
       2015 VBT Relative Risk Table - Female, 50% Non-Smoker, Age Last Birthday, Select
```

"""
function MortalityTable(lines::CSV.File)
	#what lines the table starts at
	table_starts = findall(line -> ~ismissing(line[1]) && line[1] == "Row\\Column",lines) .+ 1

	# Construct MetaData
	raw_meta = Dict()
	for line in lines
		if ismissing(line[1])
			break
		end
		raw_meta[line[1]] = line[2]
	end
	
	d = TableMetaData(
		name = get(raw_meta,"Table Name:",nothing),
		id = get(raw_meta,"Table Identity:",nothing),
		provider = get(raw_meta,"Provider Name:",nothing),
		reference = get(raw_meta,"Table Reference:",nothing),
		content_type = get(raw_meta,"Content Type:",nothing),
		description = get(raw_meta,"Table Description:",nothing),
		comments = get(raw_meta,"Comments:",nothing),
	)

# 	scale = get(raw_meta,"Scaling Factor:",nothing)
	
	# Extract values

	# figure out where table ends
	table_ends= [last_values_line(lines,ts) for ts in table_starts]


	# Parse into table
	
	if length(table_starts) == 1 
		# ultimate only
		ult_start, ult_end = table_starts[1],table_ends[1]
		ult_rates = map(lines[ult_start:ult_end]) do row
			parse(Float64,row[2])
		end
		
		ult = UltimateMortality(ult_rates,start_age=parse(Int,lines[ult_start][1]))
		
		return MortalityTable(ult; metadata=d)
	else 
		# select and ultimate
		ult_start, ult_end = table_starts[2],table_ends[2]
		ult_rates = map(lines[ult_start:ult_end]) do row
			parse(Float64,row[2])
		end
		
		ult = UltimateMortality(ult_rates,start_age=parse(Int,lines[ult_start][1]))

		sel_start, sel_end = table_starts[1],table_ends[1]

		sel_start_age = parsemaybe(Int,lines[sel_start][1])
		sel_rates = OffsetArray(
			map(sel_start:sel_end) do r
				start_age = sel_start_age + r - sel_start
				row = lines[r]
				select_rates = _select_rates([row[c] for c in 2:length(row)])
				return MortalityTables._select_row(start_age, select_rates, ult)
			end,
			sel_start_age - 1
		)
		
		return MortalityTable(sel_rates,ult,metadata=d)

	end


end

function last_values_line(lines,startline)
	for i in startline:lastindex(lines)
		if ismissing(lines[i][1]) | startswith(lines[i][1], "Table")
			return i - 1
		end
	end
	return lastindex(lines)
end

# because of the poor standardization of the CSV formatted tables from mort.SOA.org,
# sometimes the value comes through as a string, sometimes as a number when CSV.jl parses it
parsemaybe(t,x) = typeof(x) <: AbstractString ? parse(t,x) : x

# The select rates in one CSV row, given the cells after the age column: trailing
# blanks are trimmed, but a leading or interior blank is kept as `missing` so a
# row with a gap is not silently truncated.
function _select_rates(cells)
	last = findlast(!ismissing, cells)
	last === nothing && return Float64[]
	return [parsemaybe(Float64, c) for c in cells[1:last]]
end

end # module
