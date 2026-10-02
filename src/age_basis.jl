function _check_age_basis_rates(rates)
    isempty(rates) && throw(ArgumentError("mortality rates cannot be empty"))
    all(q -> zero(q) <= q <= one(q), rates) ||
        throw(ArgumentError("mortality rates must be between zero and one"))
    isone(last(rates)) ||
        throw(ArgumentError("the terminal mortality rate must be one"))
end

function _with_same_ages(rates, source)
    return UltimateMortality(rates; start_age = firstindex(source))
end

# Formula from "1980 CSO and 1980 CET Mortality Tables on an Age Last Birthday
# Basis", Transactions of the Society of Actuaries, Vol. XXXIII (1981), p. 673.

"""
    age_nearest_to_age_last(rates)

Convert age-nearest-birthday mortality rates to age-last-birthday rates,
assuming a uniform distribution of deaths. The returned vector has the same
age indices as `rates`.

The terminal rate must be one because the conversion at the final age would
otherwise require the next age's rate.
"""
function age_nearest_to_age_last(rates::AbstractVector{<:Real})
    _check_age_basis_rates(rates)
    nearest = float.(collect(rates))
    last_birthday = similar(nearest)

    for i in firstindex(nearest):(lastindex(nearest)-1)
        q = nearest[i]
        q_next = nearest[i+1]
        last_birthday[i] = (q + (one(q) - q) * q_next) / (one(q) + one(q) - q)
    end
    last_birthday[end] = nearest[end]

    return _with_same_ages(last_birthday, rates)
end

"""
    age_last_to_age_nearest(rates)

Convert age-last-birthday mortality rates to age-nearest-birthday rates,
assuming a uniform distribution of deaths. The returned vector has the same
age indices as `rates`.

The inverse is calculated backward from a terminal rate of one.
"""
function age_last_to_age_nearest(rates::AbstractVector{<:Real})
    _check_age_basis_rates(rates)
    last_birthday = float.(collect(rates))
    nearest = similar(last_birthday)
    nearest[end] = last_birthday[end]

    for i in (lastindex(last_birthday)-1):-1:firstindex(last_birthday)
        q = last_birthday[i]
        q_next = nearest[i+1]
        nearest[i] = ((one(q) + one(q)) * q - q_next) / (one(q) + q - q_next)
    end

    all(q -> zero(q) <= q <= one(q), nearest) ||
        throw(ArgumentError("rates do not produce valid age-nearest-birthday probabilities"))

    return _with_same_ages(nearest, rates)
end
