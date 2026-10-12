# Numerical kernels shared by the parametric laws (`parameterized_models.jl`, `ratio_laws.jl`).

# expm1(x)/x and log1p(x)/x, both 1 at x = 0: smooth there, so the closed forms below keep
# their zero-growth limits (and derivatives) instead of evaluating 0/0. They are used only
# near b·age = 0: at large |b·age| and infinite ages the factored forms would evaluate
# Inf/Inf or 0·Inf, so the plain closed forms apply there.
_exprel(x) = abs(x) < 1e-5 ? 1 + x / 2 + x^2 / 6 + x^3 / 24 : expm1(x) / x
_log1pdivx(x) = abs(x) < 1e-5 ? 1 - x / 2 + x^2 / 3 - x^3 / 4 : log1p(x) / x

# k·age, which is 0 for k = 0 even at an infinite age: a law's zero coefficient contributes
# nothing there (b = 0 is a constant hazard, c = 0 the Gompertz law).
_times_age(k, age) = iszero(k) ? zero(k * age) : k * age

# a / (e + a·c) for e > 0 and c ≥ 0, the form the ratio laws take once divided through by their
# exponential. Where a·c ≤ e it is evaluated as written: smooth in a at a = 0, and its derivative
# in a, e / (e + a·c)², is formed without cancellation. Where a·c > e it is also divided by a,
# 1 / (e/a + c), which stays finite where a·c overflows and keeps its derivatives free of
# cancellation and of the overflowing (e + a·c)². The two forms are the same function, so at
# a·c = e either gives the value and derivatives (a ForwardDiff comparison can break a tie between
# equal primal values by their partials, which then doesn't matter). Under ForwardDiff 1.x `iszero`
# in the laws' zero-coefficient shortcuts also requires zero partials, so a derivative at a = 0
# goes through this form.
_amplitude_ratio(a, c, e) = a * c <= e ? a / (e + a * c) : inv(e / a + c)

# a·exp(x) / (1 + k·a·exp(x)), the bounded ratio in Beard's and Kannisto's laws. The algebra
# follows the sign of x: for x > 0 it is divided through by exp(x), a / (exp(-x) + k·a), finite
# where exp(x) overflows and tending to 1/k; for x ≤ 0 the numerator z = a·exp(x) is formed
# directly, z / (1 + k·z), since there it is exp(-x) that can overflow.
function _logistic_ratio(a, k, x)
    iszero(a) && return zero(a * x)   # no hazard, even at an infinite age
    x > 0 && return _amplitude_ratio(a, k, exp(-x))
    z = a * exp(x)
    return _amplitude_ratio(z, k, one(z))
end

# log(1 + exp(y)), without overflow for large y
_log1pexp(y) = y > 0 ? y + log1p(exp(-y)) : log1p(exp(y))
