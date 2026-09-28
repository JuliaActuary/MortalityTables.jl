module MortalityTablesCSVExt

using MortalityTables, CSV
import MortalityTables: MortalityTable, TableMetaData

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
	
	# Parse each block of rates, from its first line to where it ends, into the labeled rates XTbML
	# parses to (the ages in the first column and the durations in the header row), which
	# `_table_from_labels` places, so grouped ages or blank cells keep every rate at its own age.
	source = "CSV table $(something(d.name, "(unnamed)"))"
	blocks = [(ts, last_values_line(lines, ts)) for ts in table_starts]
	select_block, ultimate_block = MortalityTables._select_and_ultimate(blocks, "tables", source)
	ultimate = _ultimate_records(lines, ultimate_block...)
	select = isnothing(select_block) ? nothing : _select_rows(lines, select_block...)
	return MortalityTables._table_from_labels((; select, ultimate, metadata = d), source)
end

# The `(age, rate)` records of the ultimate rates between two lines.
_ultimate_records(lines, first_line, last_line) =
	[(age = parsemaybe(Int, lines[r][1]), rate = parsemaybe(Float64, lines[r][2])) for r in first_line:last_line]

# The `(issue_age, rates)` records of the select rates between two lines, whose durations label the
# header row above them. A row without any rate has no select period: like an issue age absent from
# the table (as in XTbML) it is left `missing`, rather than read as a zero-length select period that
# falls straight through to the ultimate rates.
function _select_rows(lines, first_line, last_line)
	header = lines[first_line - 1]
	durations = [ismissing(header[c]) ? missing : parsemaybe(Int, header[c]) for c in 2:length(header)]
	rate_rows = [r for r in first_line:last_line if any(c -> !ismissing(lines[r][c]), 2:length(lines[r]))]
	return map(rate_rows) do r
		row = lines[r]
		(issue_age = parsemaybe(Int, row[1]), rates = _select_records([row[c] for c in 2:length(row)], durations))
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
# (a numeric column makes the header's duration labels floats, such as `4.0`)
parsemaybe(t,x) = typeof(x) <: AbstractString ? parse(t,x) : convert(t,x)

# The `(duration, rate)` records of one CSV select row, given the cells after the age column and
# the duration labels of the table's header row. A blank cell is a `missing` rate, which
# `_table_from_labels` treats as no rate: trailing blanks shorten the row and a leading or interior
# blank leaves `missing` at its duration, so a row with a gap is neither truncated nor shifted.
_select_records(cells, durations) =
	[(duration = durations[i], rate = ismissing(cells[i]) ? missing : parsemaybe(Float64, cells[i])) for i in eachindex(cells)]

end # module
