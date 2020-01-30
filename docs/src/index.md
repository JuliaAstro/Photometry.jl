```@meta
CurrentModule = Photometry
DocTestSetup = :(using Photometry)
```

# Photometry

[![Build Status](https://github.com/mileslucas/AperturePhotometry.jl/workflows/CI/badge.svg)](https://github.com/mileslucas/AperturePhotometry.jl/actions)
[![Coverage](https://codecov.io/gh/mileslucas/AperturePhotometry.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mileslucas/AperturePhotometry.jl)

## Installation

To install from the REPL, enter Pkg-mode (`]`)

```julia-repl
pkg> add https://github.com/mileslucas/Photometry.jl
```

## Getting Started

Here is a basic example to do some aperture photometry using [`CircularAperture`](@ref) and [`CircularAnnulus`](@ref).

```jldoctest
data = ones(100, 100)

ap1 = CircularAperture(50, 50, 3)
ap2 = CircularAperture(28, 32, 5)

results = aperture_photometry([ap1, ap2], data)
@assert results.aperture_sum[1] ≈ 9π
@assert results.aperture_sum[2] ≈ 25π

results

# output
2×3 DataFrames.DataFrame
│ Row │ xcenter │ ycenter │ aperture_sum │
│     │ Int64   │ Int64   │ Float64      │
├─────┼─────────┼─────────┼──────────────┤
│ 1   │ 50      │ 50      │ 28.2743      │
│ 2   │ 28      │ 32      │ 78.5398      │
```

## Contributing

If you're interested in contributing, go ahead and check out the [issues](https://github.com/mileslucas/Photometry.jl/issues) or make a [pull request](https://github.com/mileslucas/Photometry.jl/pulls). If you add a new feature, please write appropriate unit tests for it and bump the package's minor version.
