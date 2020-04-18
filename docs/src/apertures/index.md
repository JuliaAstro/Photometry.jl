# Aperture Photometry

## Introduction

Aperture photometry uses [Apertures](@ref) to cut out and sum values in an image. A very basic mask might be a square of pixels at a certain position. We can model this as a matrix of ones and zeros like
```julia
[0 0 0 0 0
 0 1 1 1 0
 0 1 1 1 0
 0 1 1 1 0
 0 0 0 0 0]
```
If we have some data like
```julia
[7 9 6 0 8
 8 5 8 7 9
 5 6 2 2 7
 9 7 3 4 1
 7 8 0 9 8]
```
then the result of our aperture photometry looks like
```julia
[0 0 0 0 0     [7 9 6 0 8     [0 0 0 0 0
 0 1 1 1 0      8 5 8 7 9      0 5 8 7 0
 0 1 1 1 0  .*  5 6 2 2 7  =   0 6 2 2 0
 0 1 1 1 0      9 7 3 4 1      0 7 3 4 0
 0 0 0 0 0]     7 8 0 9 8]     0 0 0 0 0]

sum(result) = 44
```

This module uses the above principal with common aperture shapes in a fast and precise manner, including exact overlaps between apertures and pixels.

The majority of the lifting is done with the [`aperture_photometry`](@ref) function with common shapes being described in [Apertures](@ref). It is possible to create a custom aperture by sub-typing the [`Aperture.AbstractAperture`](@ref) class, although it may be easier to perform PSF photometry instead.

## Pixel Convention

`Photometry.jl` follows the same convention as FITS, WCS, IRAF, ds9, and SourceExtractorBackground with `(1, 1)` being the _center_ on the bottom-left pixel. This means the exact bottom-left corner is at `(0.5, 0.5)`. Pixels increase up and to the right until `axis_length + 0.5`.

This is mostly in line with Julia's indexing, although it is important to remember that arrays are layed out in `(y, x)` due to the row-column interface. So the pixel at `(34, 56)` is at `image[56, end-34]`.

## API/Reference

```@docs
aperture_photometry
```

## Performance

Below is a benchmark result comparing Photometry.jl to [photutils](https://github.com/astropy/photutils). The benchmark code can be found in the [benchmarks folder](https://github.com/JuliaAstro/Photometry.jl/blob/master/benchmarks/circle).

![](circle_apertures_benchmark.png)

```julia
julia> versioninfo()
Julia Version 1.4.0
Commit b8e9a9ecc6 (2020-03-21 16:36 UTC)
Platform Info:
  OS: macOS (x86_64-apple-darwin18.6.0)
  CPU: Intel(R) Core(TM) i5-8259U CPU @ 2.30GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-8.0.1 (ORCJIT, skylake)
```
