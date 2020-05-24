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

# Reference: Experience Study Calculations, 2016, Society of Actuaries
# https://www.soa.org/globalassets/assets/Files/Research/2016-10-experience-study-calculations.pdf

"""
    q(MortalityTable,issue_age,duration,time,DeathDistribution)

A function to calculate non-whole year force of mortality, where you must
    make an assumption (`DeathDistribution`) about how the annual rate of
    mortaility applies through the rest of the year.
"""
function q(tbl, issue_age, duration, time, dist)
    return 1 - p(tbl, issue_age, duration, time, dist)
end

"""
    p(MortalityTable,issue_age,duration,time,DeathDistribution)

A function to calculate non-whole year survivorship, where you must
    make an assumption (`DeathDistribution`) about how the annual rate of
    mortaility applies through the rest of the year.
"""
function p(tbl, issue_age, duration, time, dist::Balducci)

    if time > 1
        whole_time = trunc(Int, time)
        p′ = p(tbl, issue_age, duration, whole_time)
        return p′ *
               p(tbl, issue_age, duration + whole_time, time - whole_time, dist)
    else
        # when the original time is a whole number, the remainder of time
        # can be zero, but the duration was incremented past the end of the table
        # as a result
        duration = time == 0 ? duration - 1 : duration
        q′ = q(tbl, issue_age, duration)

        return (1 - q′) / (1 - (1 - time) * q′)
    end
end

"""
    p(MortalityTable,issue_age,duration,time,DeathDistribution)

A function to calculate non-whole year survivorship, where you must
    make an assumption (`DeathDistribution`) about how the annual rate of
    mortaility applies through the rest of the year.
"""
function p(tbl, issue_age, duration, time, dist::Constant)

    if time > 1
        whole_time = trunc(Int, time)
        p′ = p(tbl, issue_age, duration, whole_time)
        return p′ *
               p(tbl, issue_age, duration + whole_time, time - whole_time, dist)
    else
        # when the original time is a whole number, the remainder of time
        # can be zero, but the duration was incremented past the end of the table
        # as a result
        duration = time == 0 ? duration - 1 : duration
        q′ = q(tbl, issue_age, duration)

        return (1 - q′)^time
    end
end

"""
    p(MortalityTable,issue_age,duration,time,DeathDistribution)

A function to calculate non-whole year survivorship, where you must
    make an assumption (`DeathDistribution`) about how the annual rate of
    mortaility applies through the rest of the year.
"""
function p(tbl, issue_age, duration, time, dist::Uniform)

    if time > 1
        whole_time = trunc(Int, time)
        p′ = p(tbl, issue_age, duration, whole_time)
        return p′ *
               p(tbl, issue_age, duration + whole_time, time - whole_time, dist)
    else
        # when the original time is a whole number, the remainder of time
        # can be zero, but the duration was incremented past the end of the table
        # as a result
        duration = time == 0 ? duration - 1 : duration
        q′ = q(tbl, issue_age, duration)

        return (1 - time * q′)
    end
end
