using Photometry
using Random
using Statistics
using ParallelTestRunner

Random.seed!(8462852)

runtests(Photometry, ARGS)
