using Documenter, MortalityTables

makedocs(;
    modules=[MortalityTables],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/JuliaActuary/MortalityTables.jl/blob/{commit}{path}#L{line}",
    sitename="MortalityTables.jl",
    authors="Alec Loudenback",
)

deploydocs(;
    repo="github.com/JuliaActuary/MortalityTables.jl",
)
