using Photometry
using Documenter

setup = quote
    using Photometry
    using Plots
    unicodeplots()
end

DocMeta.setdocmeta!(Photometry, :DocTestSetup, setup; recursive = true)

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
            "Getting Started" => "apertures/index.md",
            "Apertures" => "apertures/apertures.md",
            "Examples" => "apertures/examples.md"
        ]
    ],)

deploydocs(;
    repo = "github.com/JuliaAstro/Photometry.jl",)
