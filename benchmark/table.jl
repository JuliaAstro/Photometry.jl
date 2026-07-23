using Pkg: PackageSpec
using AirspeedVelocity
using AirspeedVelocity.TableUtils: format_time, format_memory
using PrettyTables

const PKG = "Photometry"
const REV = "dirty"
const RESULTS_DIR = joinpath(@__DIR__, "results")
const RESULTS_FILE = joinpath(RESULTS_DIR, "results_$(PKG)@$(REV).json")

if !isfile(RESULTS_FILE)
    @info "No results found at $RESULTS_FILE — running benchpkg"
    benchpkg(PKG;
        rev = REV,
        path = dirname(@__DIR__),
        output_dir = RESULTS_DIR,
        script = joinpath(@__DIR__, "benchmarks.jl"),
        dont_print = true,
    )
end

specs = [PackageSpec(; name = PKG, rev = REV)]
results = load_results(specs; input_dir = RESULTS_DIR)
benches = first(values(results))

names = collect(keys(benches))
if "time_to_load" in names
    deleteat!(names, findfirst(==("time_to_load"), names))
    push!(names, "time_to_load")
end

times = [format_time(benches[n]) for n in names]
mems  = [format_memory(benches[n]) for n in names]

pretty_table(hcat(names, times, mems);
    column_labels = ["benchmark", "time (median ± IQR)", "memory"],
    alignment = [:l, :r, :r],
    fit_table_in_display_vertically = false,
    fit_table_in_display_horizontally = false,
)
