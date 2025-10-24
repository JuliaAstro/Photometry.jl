using ParallelTestRunner: runtests
import Photometry

const init_code = quote
    import StatsBase: median, mean, std, mad

    const DATA_DIR = joinpath(@__DIR__, "data")
end

runtests(Photometry, Base.ARGS; init_code)
