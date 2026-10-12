"""
    DeathDistribution

An abstract type used to form an assumption of how deaths occur throughout a
    year. See `Balducci()`, `UniformDeaths()`, and `ConstantForce()` for concrete
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
    UniformDeaths()

A `DeathDistribution` type that assumes deaths are spread uniformly over the year of age
(UDD), which gives an increasing force of mortality over the year.
"""
struct UniformDeaths <: DeathDistribution end

"""
    ConstantForce()

A `DeathDistribution` type that assumes a constant force of mortality over the year of age.
"""
struct ConstantForce <: DeathDistribution end