using Photometry
using Documenter

setup = quote
    using Photometry
    using Random
    Random.seed!(123456)
end
DocMeta.setdocmeta!(Photometry, :DocTestSetup, setup; recursive = true)

makedocs(
    modules = [Photometry],
    authors = "Miles Lucas <mdlucas@hawaii.edu>",
    repo = "https://github.com/JuliaAstro/Photometry.jl/blob/{commit}{path}#L{line}",
    sitename = "Photometry.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://juliaastro.github.io/Photometry.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Background Estimation" => [
            "Getting Started" => "background/index.md",
            "Background Estimators" => "background/estimators.md",
            "Background Interpolators" => "background/interpolators.md"
        ],
        "Source Detection" => [
            "Getting Started" => "detection/index.md",
            "Algorithms" => "detection/algs.md",
        ],
        "Aperture Photometry" => [
            "Getting Started" => "apertures/index.md",
            "Apertures" => "apertures/apertures.md",
            "Examples" => "apertures/examples.md"
        ]
   ])

deploydocs(
    repo = "github.com/JuliaAstro/Photometry.jl",
    push_preview = true)
