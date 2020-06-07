"""
    DeathDistribution

An abstract type used to form an assumption of how deaths occur throughout a
    year. See `Balducci()`, `Uniform()`, and `Constant()` for concrete
    assumption types.
"""
abstract type DeathDistribution end

"""
    Balducci()

A `DeathDistribution` type that assumes a decreasing force of mortality
over the year.
"""
struct Balducci <: DeathDistribution end

"""
    Uniform()

A `DeathDistribution` type that assumes an increasing force of mortality
over the year.
"""
struct Uniform <: DeathDistribution end

"""
    Constant()

A `DeathDistribution` type that assumes a constant force of mortality
over the year.
"""
struct Constant <: DeathDistribution end