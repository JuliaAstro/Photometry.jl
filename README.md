# Photometry.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaastro.org/Photometry/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaastro.org/Photometry.jl/dev)

[![CI (stable)](https://github.com/JuliaAstro/Photometry.jl/actions/workflows/ci_stable.yml/badge.svg)](https://github.com/JuliaAstro/Photometry.jl/actions/workflows/ci_stable.yml)
[![CI (pre)](https://github.com/JuliaAstro/Photometry.jl/actions/workflows/ci_pre.yml/badge.svg)](https://github.com/JuliaAstro/Photometry.jl/actions/workflows/ci_pre.yml)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/P/Photometry.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![codecov](https://codecov.io/gh/JuliaAstro/Photometry.jl/graph/badge.svg?token=lqTjxxg5dg)](https://codecov.io/gh/JuliaAstro/Photometry.jl)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-orange.svg)](https://opensource.org/licenses/BSD-3-Clause)

This is a package for performing astronomical photometry using modern and efficient algorithms.

Inspired by [photutils](https://github.com/astropy/photutils), [SEP](https://github.com/kbarbary/sep), and [AperturePhotometry.jl](https://github.com/kbarbary/AperturePhotometry.jl).

## Usage

Here is a basic example to do some aperture photometry using `CircularAperture`. The `photometry` function performs the photometry using a given method. Please see [the documentation](https://JuliaAstro.github.io/Photometry.jl/dev) for more examples and reference material.

```julia
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

Output:
```plain
Table with 4 columns and 2 rows:
     xcenter  ycenter  aperture_sum  aperture_sum_err
   ┌─────────────────────────────────────────────────
 1 │ 50.0     50.0     28.2743       5.31736
 2 │ 0.5      0.5      19.635        4.43113
```

## Contributing

Please see the to-do list above for project ideas as well as any open issues! If you add new functionality, please add appropriate documentation and testing. In addition, please increment the minor version of the package to reflect the new changes!

Tests are run with [ParallelTestRunner.jl](https://github.com/JuliaTesting/ParallelTestRunner.jl). See below for some brief usage examples:

**`Pkg.test()` workflow in root directory:**

```julia-repl
julia --proj --threads=auto

(Photometry) pkg> test
```
**Interactive workflow in `test` directory:**

```julia-repl
julia --proj --threads=auto

julia> using Photometry, ParallelTestRunner
```

```julia-repl
# List available testsets

julia> runtests(Photometry, ["--list"])
Available tests:
 - aperture/circular
 - aperture/elliptical
 - aperture/overlap
 - aperture/photometry
 - aperture/plotting
 - aperture/rectangle
 - backgroud/background
 - backgroud/estimators
 - backgroud/interpolators
 - detection/detection
```

```julia-repl
# Run a subset of tests

julia> const init_code = quote
           import StatsBase: median, mean, std, mad

           const DATA_DIR = joinpath(@__DIR__, "data")
       end

julia> runtests(Photometry, ["--verbose"]; init_code, test_filter = test -> occursin("aperture", test))

Test Summary:           | Pass  Total   Time
  Overall               | 9478   9478  24.8s
    aperture/plotting   |   27     27   6.4s
    aperture/circular   |   21     21   1.6s
    aperture/overlap    | 9226   9226   9.9s
    aperture/elliptical |   20     20   1.2s
    aperture/rectangle  |   12     12   0.3s
    aperture/photometry |  172    172  19.9s
    SUCCESS
```

## License

The work derived from `astropy/photutils` is BSD 3-clause and the work derived from `kbarbary/sep` is BSD 3-clause. All other work is considered MIT expat. Therefore this work as a whole is BSD 3-clause. [`LICENSE`](LICENSE) contains all licenses and any files using derived work are noted at the top of the file.
