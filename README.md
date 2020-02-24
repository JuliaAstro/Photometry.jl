# Photometry.jl

[![CI](https://github.com/juliaastro/photometry.jl/workflows/CI/badge.svg?branch=master)](https://github.com/juliaastro/photometry.jl/workflows/CI)
[![Coverage](https://codecov.io/gh/JuliaAstro/Photometry.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaAstro/Photometry.jl)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaAstro.github.io/Photometry.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaAstro.github.io/Photometry.jl/dev)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-orange.svg)](https://opensource.org/licenses/BSD-3-Clause)

This is a package for performing astronomical photometry using modern and efficient algorithms.

Inspired by [photutils](https://github.com/astropy/photutils), [SEP](https://github.com/kbarbary/sep), and [AperturePhotometry.jl](https://github.com/kbarbary/AperturePhotometry.jl).

## To-do list

Here are features planned
- [x] Circular Aperture, Circular Annulus
- [x] Edge handling
- [x] Elliptical Aperture, Elliptical Annulus
- [ ] Rectangular Aperture, Rectangular Annulus
- [x] Plotting for aperture types
- [ ] Using SkyCoords/WCS for positions
  - Needs some work done in WCS
- [ ] background modelling
- [ ] star extraction
- [ ] filtering

In addition, the funcitonality needs documented, tested, and benchmarked. 

## Usage

Please see [the documentation](https://JuliaAstro.github.io/Photometry.jl/dev) for examples and reference material.

## Contributing

Please see the to-do list above for project ideas as well as any open issues! If you add new functionality, please add appropriate documentation and testing. In addition, please increment the minor version of the package to reflect the new changes!

## License

The work derived from `astropy/photutils` is BSD 3-clause and the work derived from `kbarbary/sep` is BSD 3-clause. All other work is considered MIT expat. Therefore this work as a whole is BSD 3-clause. [`LICENSE`](LICENSE) contains all licenses and any files using derived work are noted at the top of the file.
