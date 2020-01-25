using AperturePhotometry
using Documenter

makedocs(;
    modules=[AperturePhotometry],
    authors="Miles Lucas <mdlucas@hawaii.edu>",
    repo="https://github.com/mileslucas/AperturePhotometry.jl/blob/{commit}{path}#L{line}",
    sitename="AperturePhotometry.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mileslucas.github.io/AperturePhotometry.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mileslucas/AperturePhotometry.jl",
)
