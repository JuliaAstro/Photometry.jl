using Photometry
using Documenter

DocMeta.setdocmeta!(Photometry, :DocTestSetup, :(using Photometry); recursive = true)

makedocs(;
    modules = [Photometry],
    authors = "Miles Lucas <mdlucas@hawaii.edu>",
    repo = "https://github.com/JuliaAstro/Photometry.jl/blob/{commit}{path}#L{line}",
    sitename = "Photometry.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://juliaastro.github.io/Photometry.jl",
        assets = String[],),
    pages = [
        "Home" => "index.md",
        "Aperture Photometry" => [
            "Apertures" => "apertures.md"
        ]
    ],)

deploydocs(;
    repo = "github.com/juliaastro/Photometry.jl",)
