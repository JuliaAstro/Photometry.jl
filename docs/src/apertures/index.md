# Aperture Photometry

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

## API/Reference

```@docs
aperture_photometry
```