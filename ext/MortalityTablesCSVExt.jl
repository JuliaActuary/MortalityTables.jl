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


	# Parse into table: rates are placed by their labels (the ages in the first column and the
	# durations in the header row), as for XTbML, so grouped ages or blank cells keep every rate
	# at its own age.
	source = "CSV table $(something(d.name, "(unnamed)"))"
	
	if length(table_starts) == 1 
		# ultimate only
		ult = _ultimate(lines, table_starts[1], table_ends[1], source)
		
		return MortalityTable(ult; metadata=d)
	else 
		# select and ultimate
		ult = _ultimate(lines, table_starts[2], table_ends[2], source)

		sel_start, sel_end = table_starts[1],table_ends[1]
		header = lines[sel_start - 1]
		durations = [ismissing(header[c]) ? missing : parsemaybe(Int, header[c]) for c in 2:length(header)]
		issue_ages = [parsemaybe(Int, lines[r][1]) for r in sel_start:sel_end]
		rows = map(issue_ages, sel_start:sel_end) do issue_age, r
			row = lines[r]
			select_rates = _select_rates(
				[row[c] for c in 2:length(row)], durations,
				"the select durations for issue age $issue_age", source
			)
			return MortalityTables._select_row(issue_age, select_rates, ult)
		end
		first_issue_age, sel = MortalityTables._by_label(issue_ages, rows, "the select issue ages", source)
		
		return MortalityTable(OffsetArray(sel, first_issue_age - 1),ult,metadata=d)

	end


end

# The ultimate rates between two lines, placed by the ages in the first column.
function _ultimate(lines, first_line, last_line, source)
	ages = [parsemaybe(Int, lines[r][1]) for r in first_line:last_line]
	rates = [parsemaybe(Float64, lines[r][2]) for r in first_line:last_line]
	start_age, placed = MortalityTables._by_label(ages, rates, "the ultimate ages", source)
	return UltimateMortality(placed, start_age = start_age)
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
# (a numeric column makes the header's duration labels floats, such as `4.0`)
parsemaybe(t,x) = typeof(x) <: AbstractString ? parse(t,x) : convert(t,x)

# The select rates in one CSV row, given the cells after the age column and the duration
# labels of the table's header row. A blank cell gives no rate: trailing blanks shorten the row,
# and a leading or interior blank leaves `missing` at its duration, so a row with a gap is
# neither truncated nor shifted.
function _select_rates(cells, durations, what = "the select durations", source = "CSV table")
	given = findall(!ismissing, cells)
	isempty(given) && return Float64[]
	_, rates = MortalityTables._by_label(
		[durations[i] for i in given], [parsemaybe(Float64, cells[i]) for i in given], what, source;
		first = 1
	)
	return rates
end

end # module
