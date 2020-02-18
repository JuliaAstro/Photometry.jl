```@meta
DocTestSetup = :(using Photometry)
```

# Apertures

All apertures will rely on a position and the shape parameters.
```julia
aperture = Aperture(x0, y0, shape_params...)
```
The position can be pixels or sky coordinates. The sky coordinate positions utilize [SkyCoords.jl](https://juliaastro.github.io/SkyCoords.jl/stable) and [WCS.jl](https://juliaastro.github.io/WCS.jl/stable) for conversion.

!!! note
    The pixel positions for these apertures follow traditional image position with 1-based indexing. This means the origin is at top-left and has index `(0.5, 0.5)` at the top-left corner and `(1, 1)` at the center.

## Circular Apertures

These apertures are parametrized by radius.

```@docs
CircularAperture
CircularAnnulus
```

## Elliptical Apertures

These apertures are parametrized by the semi-major axis `a`, semi-minor axis `b` and rotation angle in degrees counter-clockwise from the positive x-axis `theta`


```@docs
EllipticalAperture
```


## Rectangular Apertures

These apertures are parametrized by side-length `a` and side-length `b`.


## API/Reference

```@docs
Aperture.AbstractAperture
mask
cutout
```
