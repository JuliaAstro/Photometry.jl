using Photometry
using Documenter
using Documenter.Remotes: GitHub

setup = quote
    using Photometry
    using Random
    Random.seed!(123456)
end

DocMeta.setdocmeta!(Photometry, :DocTestSetup, setup; recursive = true)

include("pages.jl")

makedocs(;
    modules = [Photometry],
    authors = "Miles Lucas <mdlucas@hawaii.edu>",
    repo = GitHub("JuliaAstro/Photometry.jl"),
    sitename = "Photometry.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://juliaastro.org/Photometry/stable/",
    ),
    pages = pages,
)

deploydocs(;
    repo = "github.com/JuliaAstro/Photometry.jl",
    push_preview = true,
    versions = ["stable" => "v^", "v#.#"], # Restrict to minor releases
)
