using ParallelTestRunner: runtests, find_tests, parse_args
import Photometry

const init_code = quote
    import StatsBase: median, mean, std, mad

    const DATA_DIR = joinpath(@__DIR__, "data")
end

args = parse_args(Base.ARGS)
testsuite = find_tests(@__DIR__)

runtests(Photometry, args; testsuite, init_code)
