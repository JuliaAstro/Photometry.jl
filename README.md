# Photometry.jl

[![Build Status](https://github.com/JuliaAstro/Photometry.jl/workflows/CI/badge.svg?branch=main)](https://github.com/JuliaAstro/Photometry.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/P/Photometry.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/JuliaAstro/Photometry.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaAstro/Photometry.jl)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-orange.svg)](https://opensource.org/licenses/BSD-3-Clause)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaAstro.github.io/Photometry.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaAstro.github.io/Photometry.jl/dev)

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

## License

The work derived from `astropy/photutils` is BSD 3-clause and the work derived from `kbarbary/sep` is BSD 3-clause. All other work is considered MIT expat. Therefore this work as a whole is BSD 3-clause. [`LICENSE`](LICENSE) contains all licenses and any files using derived work are noted at the top of the file.
