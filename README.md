# Photometry.jl

[![Build Status](https://github.com/JuliaAstro/Photometry.jl/workflows/CI/badge.svg)](https://github.com/JuliaAstro/Photometry.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaAstro/Photometry.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaAstro/Photometry.jl)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaAstro.github.io/Photometry.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaAstro.github.io/Photometry.jl/dev)

This is a package for performing aperture photometry and PSF photometry.

Heavily inspired by [photutils](https://github.com/astropy/photutils) and [kbarbary/AperturePhotometry.jl](https://github.com/kbarbary/AperturePhotometry.jl).

## To-do list

Here are features planned
- [x] Circular Aperture, Circular Annulus
- [x] Edge handling
- [ ] Elliptical Aperture, Elliptical Annulus
- [ ] Rectangular Aperture, Rectangular Annulus
- [x] Plotting for aperture types
- [ ] Using SkyCoords/WCS for positions
  - Needs some work done in WCS

In addition, the funcitonality needs documented, tested, and benchmarked. 

## Usage

Please see [the documentation](https://JuliaAstro.github.io/Photometry.jl/dev) for examples and reference material.

## Contributing

Please see the to-do list above for project ideas as well as any open issues! If you add new functionality, please add appropriate documentation and testing. In addition, please increment the minor version of the package to reflect the new changes!
