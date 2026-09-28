# Table files label every rate with its age or duration (XTbML elements, CSV row and column
# headers). Each value goes to its label, so a table whose labels skip (rates at grouped ages,
# or an empty cell inside a select row) keeps its values at the right ages, with `missing` where
# it gives no rate. The element type widens only when there are such gaps. A repeated label is
# ambiguous and throws; `source` names the table in that error. `first` is the label of the
# first position (durations start at 1; ages at the smallest label).
function _by_label(labels, values, what, source; first = minimum(labels))
    allunique(labels) || throw(
        ArgumentError("$source: $what repeat a label, so their rates cannot be placed: $(labels)")
    )
    n = maximum(labels) - first + 1
    labels == first:(first + n - 1) && return first, collect(values)
    placed = Vector{Union{Missing, eltype(values)}}(missing, n)
    for (label, value) in zip(labels, values)
        placed[label - first + 1] = value
    end
    return first, placed
end

# A table from its labeled rates, as the file readers (XTbML, and CSV through its extension) parse
# them: `ultimate` is a vector of `(age, rate)`; `select` is `nothing` for an ultimate-only table, or
# a vector of `(issue_age, rates)` where `rates` is a vector of `(duration, rate)` for the durations
# the file gives. Every rate goes to its label (`_by_label`), each select row takes the ultimate
# rates after its select period (`_select_row`), and `source` names the table in errors.
function _table_from_labels(tbl, source)
    start_age, ult_rates = _by_label([v.age for v in tbl.ultimate], [v.rate for v in tbl.ultimate], "the ultimate ages", source)
    ult = UltimateMortality(ult_rates, start_age = start_age)
    isnothing(tbl.select) && return MortalityTable(ult, metadata = tbl.metadata)
    rows = map(tbl.select) do (issue_age, rates)
        _, select_rates = _by_label(
            [r.duration for r in rates], [r.rate for r in rates],
            "the select durations for issue age $issue_age", source; first = 1
        )
        return _select_row(issue_age, select_rates, ult)
    end
    first_issue_age, sel = _by_label([r.issue_age for r in tbl.select], rows, "the select issue ages", source)
    return MortalityTable(OffsetArray(sel, first_issue_age - 1), ult, metadata = tbl.metadata)
end
