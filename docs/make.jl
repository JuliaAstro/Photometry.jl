using Photometry
using Documenter

setup = quote
    using Photometry
    using Random
    Random.seed!(123456)
end
DocMeta.setdocmeta!(Photometry, :DocTestSetup, setup; recursive = true)

include("pages.jl")
makedocs(
    modules = [Photometry],
    authors = "Miles Lucas <mdlucas@hawaii.edu>",
    repo = "https://github.com/JuliaAstro/Photometry.jl/blob/{commit}{path}#L{line}",
    sitename = "Photometry.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://juliaastro.github.io/Photometry.jl",
    ),
    pages=pages
)

deploydocs(
    repo = "github.com/JuliaAstro/Photometry.jl",
    push_preview = true)
