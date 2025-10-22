# Photometry.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaastro.org/Photometry/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaastro.org/Photometry.jl/dev)

[![CI](https://github.com/JuliaAstro/Photometry.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaAstro/Photometry.jl/actions/workflows/CI.yml)
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

Tests are run with [ParallelTestRunner.jl](https://github.com/JuliaTesting/ParallelTestRunner.jl), which hooks into the standard `Pkg.test` harness for more granular control. See below for some brief usage examples run from the package root directory:

```julia-repl
julia --proj

julia> using Pkg
```

**List available testsets:**

```julia-repl
julia> Pkg.test("Photometry"; test_args=`--list`);
# ...
     Testing Running tests...
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
     Testing Photometry tests passed
```

**Run a subset of testsets in verbose mode, and enable default threading:**

```julia-repl
julia> Pkg.test("Photometry"; test_args=`--verbose aperture`, julia_args=`--threads=auto`);
# ...
     Testing Running tests...
Running 3 tests in parallel. If this is too many, specify the `--jobs=N` argument to the tests, or set the `JULIA_CPU_THREADS` environment variable.
                          │          │ ──────────────── CPU ──────────────── │
Test             (Worker) │ Time (s) │ GC (s) │ GC % │ Alloc (MB) │ RSS (MB) │
aperture/photometry   (1) │        started at 2025-10-21T16:36:40.959
aperture/overlap      (2) │        started at 2025-10-21T16:36:40.959
aperture/plotting     (3) │        started at 2025-10-21T16:36:40.959
aperture/plotting     (3) │     6.33 │   0.59 │  9.3 │     638.40 │   733.98 │
aperture/circular     (3) │        started at 2025-10-21T16:36:49.971
aperture/overlap      (2) │     8.83 │   0.61 │  6.9 │     744.02 │   733.98 │
aperture/elliptical   (2) │        started at 2025-10-21T16:36:52.293
aperture/circular     (3) │     1.61 │   0.00 │  0.0 │     160.19 │   733.98 │
aperture/rectangle    (3) │        started at 2025-10-21T16:36:52.571
aperture/elliptical   (2) │     0.74 │   0.00 │  0.0 │      79.68 │   733.98 │
aperture/rectangle    (3) │     1.56 │   0.00 │  0.0 │     142.79 │   734.75 │
aperture/photometry   (1) │    19.37 │   0.94 │  4.8 │    2712.83 │   821.66 │

Test Summary:           | Pass  Total   Time
  Overall               | 9478   9478  24.1s
    aperture/plotting   |   27     27   6.3s
    aperture/overlap    | 9226   9226   8.8s
    aperture/circular   |   21     21   1.6s
    aperture/elliptical |   20     20   0.7s
    aperture/rectangle  |   12     12   1.6s
    aperture/photometry |  172    172  19.4s
    SUCCESS
     Testing Photometry tests passed
```

## License

The work derived from `astropy/photutils` is BSD 3-clause and the work derived from `kbarbary/sep` is BSD 3-clause. All other work is considered MIT expat. Therefore this work as a whole is BSD 3-clause. [`LICENSE`](LICENSE) contains all licenses and any files using derived work are noted at the top of the file.
