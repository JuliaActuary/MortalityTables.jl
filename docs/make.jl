using Documenter, MortalityTables

makedocs(;
    modules=[MortalityTables],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Usage and Examples" => "examples.md",
            "Tables" => "Tables.md",
            "Parametric Models" => "ParametricMortalityModels.md",
            "Comparison Tool" => "ComparisonTool.md",
            "Package Design" => "design.md",
        ],
        "Reference" => "Reference.md",
    ],
    repo="https://github.com/JuliaActuary/MortalityTables.jl/blob/{commit}{path}#L{line}",
    sitename="MortalityTables.jl",
    authors="Alec Loudenback",
)

deploydocs(;
    repo="github.com/JuliaActuary/MortalityTables.jl",
)
