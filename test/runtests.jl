using ParallelTestRunner: runtests, parse_args
import Photometry

const init_code = quote
    import StatsBase: median, mean, std, mad

    const DATA_DIR = joinpath(@__DIR__, "data")
end

args = parse_args(Base.ARGS)

runtests(Photometry, args; init_code)
