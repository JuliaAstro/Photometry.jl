```@meta
CurrentModule = Photometry
DocTestSetup = :(using Photometry)
```

# Photometry.jl

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/juliaastro/Photometry.jl)
[![CI](https://github.com/JuliaAstro/Photometry.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaAstro/Photometry.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaAstro/Photometry.jl/graph/badge.svg?token=lqTjxxg5dg)](https://codecov.io/gh/JuliaAstro/Photometry.jl)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-orange.svg)](https://opensource.org/licenses/BSD-3-Clause)

## Installation

To install from the REPL, enter Pkg-mode (`]`)

```julia-repl
pkg> add Photometry
```

## Getting Started

Here is a basic example to do some aperture photometry using [`CircularAperture`](@ref). The [`photometry`](@ref) function performs the photometry using a given method.

```@example
using Photometry
data = ones(100, 100)
err = ones(100, 100)

ap1 = CircularAperture(50, 50, 3)
# partial overlap
ap2 = CircularAperture(0.5, 0.5, 5)

results = photometry([ap1, ap2], data, err)
@assert results.aperture_sum[1] ≈ 9π
@assert results.aperture_sum[2] ≈ 25π / 4

results
```

## Contributing

If you're interested in contributing, go ahead and check out the [issues](https://github.com/juliaastro/Photometry.jl/issues) or make a [pull request](https://github.com/juliaastro/Photometry.jl/pulls). If you add a new feature, please write appropriate unit tests for it and bump the package's minor version.

## License

The work derived from `astropy/photutils` is BSD 3-clause and the work derived from `kbarbary/sep` is BSD 3-clause. All other work is considered MIT expat. Therefore this work as a whole is BSD 3-clause. [`LICENSE`](https://github.com/JuliaAstro/Photometry.jl/blob/main/LICENSE) contains all licenses and any files using derived work are noted at the top of the file.
