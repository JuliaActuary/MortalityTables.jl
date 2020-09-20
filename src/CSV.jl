""" 
    readcsv(path)

Read and parse a CSV file in the SOA's usual CSV format. 

# Examples

```julia-repl
julia> readcsv("path/to/table.csv")
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
function readcsv(path)
	lines = CSV.File(path,silencewarnings=true,select=[1,2],header=false)
	two_cols = [(x[1],x[2]) for x in lines]



	#what lines the table starts at
	table_starts = findall(line -> ~ismissing(line[1]) && line[1] == "Row\\Column",two_cols) .+ 1

	# Construct MetaData

	raw_meta = Dict(lines[row][1] => lines[row][2] for row in 1:table_starts[1]-1)

	
	d = TableMetaData(
		name = get(raw_meta,"Table Name:",nothing),
		id = get(raw_meta,"Table Identity:",nothing),
		provider = get(raw_meta,"Provider Name:",nothing),
		reference = get(raw_meta,"Table Reference:",nothing),
		content_type = get(raw_meta,"Content Type:",nothing),
		description = get(raw_meta,"Table Description:",nothing),
		comments = get(raw_meta,"Comments:",nothing),
		source_path = path,
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
			Parsers.parse(Float64,row[2])
		end
		
		ult = UltimateMortality(ult_rates,start_age=Parsers.parse(Int,lines[ult_start][1]))
		
		return MortalityTable(ult,d)
	else 
		# select and ultimate
		ult_start, ult_end = table_starts[2],table_ends[2]
		ult_rates = map(lines[ult_start:ult_end]) do row
			Parsers.parse(Float64,row[2])
		end
		
		ult = UltimateMortality(ult_rates,start_age=Parsers.parse(Int,lines[ult_start][1]))
		
		
		sel_start, sel_end = table_starts[1],table_ends[1]
		sel_rates 	= [
			ismissing(lines[r][c]) ? missing : Parsers.parse(Float64,lines[r][c]) for r in sel_start:sel_end, c in 2:length(lines[1])
			
			]
		sel = SelectMortality(sel_rates,ult,start_age=Parsers.parse(Int,lines[sel_start][1]))
		
		return MortalityTable(sel,ult,metadata=d)

	end


end